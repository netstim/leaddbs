function [fibsvalBin, fibsvalSum, fibsvalMean, fibsvalPeak, fibsval5Peak, fibcell, connFiberInd, totalFibers] = ea_discfibers_calcvals(vatlist, cfile, thresh)
% Calculate fiber connection values based on the VATs and the connectome

disp('Load Connectome...');
load(cfile, 'fibers', 'idx');

prefs = ea_prefs;
if ~exist('thresh','var')
    thresh = prefs.machine.vatsettings.horn_ethresh*1000;
end
[numPatient, numSide] = size(vatlist);

fibsvalBin = cell(1, numSide);
fibsvalSum = cell(1, numSide);
fibsvalMean = cell(1, numSide);
fibsvalPeak = cell(1, numSide);
fibsval5Peak = cell(1, numSide);

fibcell = cell(1, numSide);
connFiberInd = cell(1, numSide);

totalFibers = length(idx); % total number of fibers in the connectome to work with global indices

for side = 1:numSide
    fibsvalBin{side} = zeros(length(idx), numPatient);
    fibsvalSum{side} = zeros(length(idx), numPatient);
    fibsvalMean{side} = zeros(length(idx), numPatient);
    fibsvalPeak{side} = zeros(length(idx), numPatient);
    fibsval5Peak{side} = zeros(length(idx), numPatient);

    disp(['Calculate for side ', num2str(side), ':']);
    for pt = 1:numPatient
        disp(['VAT ', num2str(pt, ['%0',num2str(numel(num2str(numPatient))),'d']), '/', num2str(numPatient), '...']);
        if isstruct(vatlist) % direct nifti structs supplied
            vat = vatlist(pt,side);
        elseif iscell(vatlist) % filenames
            if strcmp(vatlist{pt,side},"skip")
                % no stimulation for this hemisphere
                continue
            elseif isfile(vatlist{pt,side})
                vat = ea_load_nii(vatlist{pt,side});
            else
                ea_cprintf('CmdWinWarnings', 'Skipping calculating connectivity: VTA doesn''t exist!\n');
                continue;
            end
        end
        % Threshold the vat efield
        vatInd = find(abs(vat.img(:))>thresh);

        % Trim connectome fibers
        [xvox, yvox, zvox] = ind2sub(size(vat.img), vatInd);
        vatmm = ea_vox2mm([xvox, yvox, zvox], vat.mat);
        filter = all(fibers(:,1:3)>=min(vatmm),2) & all(fibers(:,1:3)<=max(vatmm), 2);

        % Skip further calculation in case VAT is totally not connected
        if ~any(filter)
            continue;
        end

        trimmedFiber = fibers(filter,:);

        % Map mm connectome fibers into VAT voxel space
        [trimmedFiberInd, ~, trimmedFiberID] = unique(trimmedFiber(:,4), 'stable');

        % older less performant code line:
        %       tic  fibVoxInd = splitapply(@(fib) {ea_mm2uniqueVoxInd(fib, vat)}, trimmedFiber(:,1:3), trimmedFiberID); 

        % replaced with new more performant code:
        % --- Vectorized mm -> voxel linear indices for ALL points ---
        

        Pmm = trimmedFiber(:,1:3);
        nP  = size(Pmm,1);

        % --- mm -> vox indices (vectorized) ---
        % Precompute once outside this function if you can:
        % MinvT = inv(vat.mat)';   % or better: MinvT = vat.mat \ eye(4); MinvT = MinvT';
        MinvT = inv(vat.mat).';    % simplest, but move it out of any loop

        V = [Pmm, ones(nP,1)] * MinvT;    % faster than "/ vat.mat.'"
        V = round(V(:,1:3));

        sz = size(vat.img);
        in = V(:,1)>=1 & V(:,1)<=sz(1) & ...
            V(:,2)>=1 & V(:,2)<=sz(2) & ...
            V(:,3)>=1 & V(:,3)<=sz(3);

        % --- drop fibers that have ANY out-of-bounds point (matches your NaN logic) ---
        % Determine per-point fiber id (must be positive integers for accumarray; if not, remap)
        fid = trimmedFiberID(:);

        % If fid is not 1..K contiguous, remap once:
        % [fid,~,fid] = unique(fid,'stable');

        badFib = accumarray(fid, ~in, [], @any, false);  % logical per fiber
        goodPoint = in & ~badFib(fid);

        % Keep only good points
        fid = fid(goodPoint);
        V   = V(goodPoint,:);

        % Also keep consistent trimming index list
        % (trimmedFiberInd is per fiber, so remove bad fibers there)
        trimmedFiberInd(badFib) = [];

        % Compute linear indices (int32 is fine; uint32 also ok)
        lin = sub2ind(sz, V(:,1), V(:,2), V(:,3));
        lin = uint32(lin);

        % --- one global sort by (fiber, lin) ---
        [~,ord] = sortrows([double(fid), double(lin)], [1 2]);  % double only for sort keys
        fid = fid(ord);
        lin = lin(ord);

        % --- remove duplicates WITHIN each fiber (replaces per-fiber unique) ---
        sameAsPrev = [false; (fid(2:end)==fid(1:end-1)) & (lin(2:end)==lin(1:end-1))];
        fid = fid(~sameAsPrev);
        lin = lin(~sameAsPrev);

        % --- split to cell array without splitapply ---
        % Find fiber boundaries in the sorted fid list
        isNew = [true; fid(2:end)~=fid(1:end-1)];
        startIdx = find(isNew);
        endIdx   = [startIdx(2:end)-1; numel(fid)];

        fibVoxInd = arrayfun(@(s,e) lin(s:e), startIdx, endIdx, 'UniformOutput', false);
        % end replacement

        % Remove outliers
        fibVoxInd(cellfun(@(x) any(isnan(x)), fibVoxInd)) = [];
        trimmedFiberInd(cellfun(@(x) any(isnan(x)), fibVoxInd)) = [];

        % Find connected fibers
        % old less performant code
        %tic; connected = cellfun(@(fib) any(ismember(fib, vatInd)), fibVoxInd); toc

        % replaced with more performant code:
        mask = false(numel(vat.img),1);
        mask(vatInd) = true;

        connected = cellfun(@(fib) any(mask(fib)), fibVoxInd);
        % end replace

        % Generate binary fibsval for the T-test method
        fibsvalBin{side}(trimmedFiberInd(connected), pt)=1;

        % Checck intersection between vat and the connected fibers
        vals = cellfun(@(fib) vat.img(intersect(fib, vatInd)), fibVoxInd(connected), 'Uni', 0);

        % Generate fibsval for the Spearman's correlation method
        fibsvalSum{side}(trimmedFiberInd(connected), pt) = cellfun(@sum, vals);
        fibsvalMean{side}(trimmedFiberInd(connected), pt) = cellfun(@mean, vals);
        fibsvalPeak{side}(trimmedFiberInd(connected), pt) = cellfun(@max, vals);
        fibsval5Peak{side}(trimmedFiberInd(connected), pt) = cellfun(@(x) mean(maxk(x,ceil(0.05*numel(x)))), vals);
    end

    % Remove values for not connected fibers, convert to sparse matrix
    fibIsConnected = any(fibsvalBin{side}, 2);
    fibsvalBin{side} = sparse(fibsvalBin{side}(fibIsConnected, :));
    fibsvalSum{side} = sparse(fibsvalSum{side}(fibIsConnected, :));
    fibsvalMean{side} = sparse(fibsvalMean{side}(fibIsConnected, :));
    fibsvalPeak{side} = sparse(fibsvalPeak{side}(fibIsConnected, :));
    fibsval5Peak{side} = sparse(fibsval5Peak{side}(fibIsConnected, :));

    % Extract connected fiber cell
    connFiberInd{side} = find(fibIsConnected);
    connFiber = fibers(ismember(fibers(:,4), connFiberInd{side}), 1:3);
    fibcell{side} = mat2cell(connFiber, idx(connFiberInd{side}));
end
