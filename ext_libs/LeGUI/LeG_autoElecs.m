function [WC,T] = LeG_autoElecs(app)

Img         = app.CTImg;
MaxElecs    = app.detect.maxelecs;

ElecRad     = 2;                 % electrode radius in mm (standard intracranial)
XYZScale    = app.XYZScale;      % voxel scaling (used earlier; now mostly for info)
MRInfo      = app.MRInfo;
ProjSurfRaw = app.ProjSurfRaw;

TMax = (app.CTRng(4)-app.CTRng(1)) / (app.CTRng(2)-app.CTRng(1)); % near max
TMin = 1; % 99th percentile
Thresh = linspace(TMax, TMin, 21);
Thresh(1) = [];

ThreshHU = Thresh*(app.CTRng(2)-app.CTRng(1)) + app.CTRng(1);
ThreshHU = (ThreshHU - app.CTInfo.pinfo(2)) ./ app.CTInfo.pinfo(1);

voxelSize = sqrt(sum(app.CTInfo.mat(1:3,1:3).^2));  % [dx, dy, dz] in mm
ElecVolVox = ceil(4/3 * pi * mean(ElecRad ./ XYZScale).^3);
minVox = max(3, round(6 * 0.4 / voxelSize(1)));

maxVox = max(ElecVolVox, minVox * 4);
ElecVolRng = [minVox, maxVox];

StartTime = tic;
CCList = cell(length(Thresh),1);
WCList = cell(length(Thresh),1);

for k = 1:length(Thresh)

    ImgBin = Img > Thresh(k);

    CC = bwconncomp(ImgBin, 6);
    CCSize = cellfun(@length, CC.PixelIdxList);

    idxRemove = CCSize < ElecVolRng(1) | CCSize > ElecVolRng(2);
    CC.PixelIdxList(idxRemove) = [];
    CC.NumObjects = numel(CC.PixelIdxList);

    if CC.NumObjects > 0
        % Compute weighted centroids and mean intensity within each CC
        s = regionprops3(CC, Img .* ImgBin, 'weightedcentroid', 'meanintensity');

        WC = s.WeightedCentroid;

        % Swap x/y to match expected coordinate conventions
        WC(:, [1,2]) = WC(:, [2,1]);

        % Convert voxel coordinates to mm using MRInfo.mat
        WCmm = [WC, ones(size(WC,1),1)] * MRInfo.mat';
        WCmm(:,4) = [];

        % Keep only detections within the projection surface (e.g., brain/skull mask)
        idxKeep = LeG_intriangulation(ProjSurfRaw.vertices, ProjSurfRaw.faces, WCmm);

        s(~idxKeep, :) = [];
        WC(~idxKeep, :) = [];
        CC.PixelIdxList(~idxKeep) = [];
        CC.NumObjects = numel(CC.PixelIdxList);

        if CC.NumObjects > 0
            % Sort detections by mean intensity (brightest first)
            [~, sidx] = sort([s.MeanIntensity], 'descend');
            WC = round(WC(sidx, :));

            WCList{k} = WC;
            CCList{k} = CC;
        else
            WCList{k} = [];
            CCList{k} = CC;
        end
    else
        WCList{k} = [];
        CCList{k} = CC;
    end
end

NumObj = cellfun(@(x) size(x,1), WCList);

if all(NumObj == 0)
    idx = round(numel(Thresh)/2);
else
    relTol    = 0.2;
    absTolMin = 2;

    allowedDelta = max(absTolMin, round(relTol * max(NumObj)));

    % Minimum number of objects considered "plausible" electrodes
    % e.g. ~30% of MaxElecs, but never less than 5
    minObjs = max(5, round(0.3 * MaxElecs));

    stableMask = abs(diff(NumObj)) <= allowedDelta & NumObj(1:end-1) >= minObjs;
    cc = bwconncomp(stableMask);

    skipflag = true;
    if cc.NumObjects > 0
        ccsize = cellfun(@length, cc.PixelIdxList);
        cc.PixelIdxList(ccsize < 2) = []; % need at least 3 points (2 diffs)
        cc.NumObjects = numel(cc.PixelIdxList);

        if cc.NumObjects > 0
            ccval = cellfun(@(x) mean(NumObj(x)), cc.PixelIdxList);
            [~, midx] = max(ccval);
            idxSeg = cc.PixelIdxList{midx};
            idx = idxSeg(round(end/2));     % pick middle threshold of that segment
            skipflag = false;
        end
    end

    % Fallback if no stable segments were found
    if skipflag
        cand = find(NumObj >= minObjs & NumObj <= MaxElecs);
        if ~isempty(cand)
            idx = cand(round(end/2));
        else
            idx = round(numel(Thresh)/2);
        end
    end
end

WC   = WCList{idx};
THU  = ThreshHU(idx);
T    = Thresh(idx);

% Outlier removal and gap filling along electrode shaft
if ~isempty(WC)
    WCmm = [WC, ones(size(WC,1),1)] * MRInfo.mat';
    WCmm(:,4) = [];

    % Remove duplicates: contacts closer than 1 mm to a neighbor
    pd0 = pdist2(WCmm, WCmm, 'euclidean', 'smallest', 2);
    nn0 = pd0(end,:);             % nearest neighbor distance (mm) for each contact
    keep = nn0 >= 1;
    WC   = WC(keep,:);
    WCmm = WCmm(keep,:);

    n = size(WCmm,1);
    if n >= 3
        shaftRadius = 8;  % tune this if shafts are closer/farther

        distMat = pdist2(WCmm, WCmm);
        adj = distMat <= shaftRadius & distMat > 0;

        if any(adj(:))
            G = graph(adj);
            clusterId = conncomp(G)';
        else
            clusterId = (1:n)';
        end
        clusterCount = max(clusterId);

        newMM = [];
        for c = 1:clusterCount
            idxC = find(clusterId == c);
            if numel(idxC) < 3
                continue;  % need at least 3 points to define spacing reliably
            end

            X  = WCmm(idxC,:);
            mu = mean(X,1);
            X0 = X - mu;

            [~,~,V] = svd(X0,'econ');
            dir = V(:,1);

            % Project onto line
            t = X0 * dir;
            [tSort, order] = sort(t);
            if numel(tSort) < 2
                continue;
            end
            gaps = diff(tSort);
            baseSpacing = median(gaps);

            if baseSpacing <= 0
                continue;
            end

            % Fill gaps that are "too big" (likely missed contacts)
            for g = 1:numel(gaps)
                gap = gaps(g);
                if gap > 1.2*baseSpacing && gap < 8*baseSpacing
                    nMissing = round(gap/baseSpacing) - 1;
                    if nMissing < 1
                        continue;
                    end

                    tStart = tSort(g);
                    tEnd   = tSort(g+1);

                    for m = 1:nMissing
                        alpha = m/(nMissing+1);
                        tNew  = tStart + alpha*(tEnd - tStart);
                        xNew  = mu + tNew * dir';

                        dToExisting = min(pdist2(xNew, WCmm));
                        if dToExisting >= 1
                            newMM = [newMM; xNew]; %#ok<AGROW>
                        end
                    end
                end
            end
        end

        if ~isempty(newMM)
            mmH  = [newMM, ones(size(newMM,1),1)]';
            voxH = MRInfo.mat \ mmH;
            newVox = voxH(1:3,:)';

            sz = size(Img);
            inMask = newVox(:,1) >= 1 & newVox(:,1) <= sz(1) & ...
                     newVox(:,2) >= 1 & newVox(:,2) <= sz(2) & ...
                     newVox(:,3) >= 1 & newVox(:,3) <= sz(3);
            newVox = newVox(inMask,:);

            newVox = round(newVox);

            WC   = [WC; newVox];
            WCmm = [WCmm; newMM(inMask,:)];

            pdAll = pdist2(WCmm, WCmm, 'euclidean', 'smallest', 2);
            nnAll = pdAll(end,:);
            keep2 = nnAll >= 1;
            WC    = WC(keep2,:);
            WCmm  = WCmm(keep2,:);
        end
    end

    if size(WC,1) > MaxElecs
        WC = WC(1:MaxElecs,:);
    end
end

% 8. Diagnostic plot: #detections vs threshold
fH = figure('position',[50,50,400,400], ...
            'name', app.PatientIDStr, ...
            'visible','off');
aH = axes('parent',fH);
plot(aH, ThreshHU, NumObj, '-');
hold(aH, 'on');
plot(aH, ThreshHU(idx), NumObj(idx), 'or');

EndTime = toc(StartTime);
xlabel(aH,'Thresh(HU)');
ylabel(aH,'#Elecs');

TMaxHU = (app.CTRng(4)-app.CTInfo.pinfo(2)) ./ app.CTInfo.pinfo(1);
TMinHU = (app.CTRng(3)-app.CTInfo.pinfo(2)) ./ app.CTInfo.pinfo(1);
title(aH, sprintf('%0.1fs, %0.0fhu (%0.0f,%0.0f)', ...
      EndTime, THU, TMinHU, TMaxHU)); %time, threshold(HU), min(HU), max(HU)

print(fH, fullfile(app.SaveDir,'AutoElec.png'),'-dpng','-r300')
close(fH);

end