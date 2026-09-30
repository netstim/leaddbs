function [WC, T] = LeG_autoElecs_countguided(app, Img, MRInfo, ProjSurfRaw, ...
    voxDim, ThrR, ThrHU, TMax, TMin, StartTime, varargin)

mode = 'chainfit';
if nargin >= 11 && ~isempty(varargin{1})
    if ischar(varargin{1}) || isstring(varargin{1})
        mode = lower(char(varargin{1}));
    end
end
expected.names = cell(0, 1);
expected.counts = zeros(0, 1);
if strcmp(mode, 'countguided')
    expected = getExpectedCountsCountGuided(app);
end
useCountGuidance = strcmp(mode, 'countguided') && ~isempty(expected.counts);
useChainOptimization = strcmp(mode, 'chainoptimization');

maxElecs = 250;
if isfield(app, 'detect') && isfield(app.detect, 'maxelecs') && ~isempty(app.detect.maxelecs)
    maxElecs = double(app.detect.maxelecs(1));
end

blob = computeBlobScoreMapCountGuided(Img, voxDim);
[candMM, candConf, candInside, candNear, candBlob, candInt] = ...
    extractCandidatesCountGuided(Img, blob, MRInfo.mat, ProjSurfRaw);
surfaceSupport = 0;
if ~isempty(candMM)
    surfaceSupport = mean(candInside | candNear);
end

fprintf('[AutoElec] Candidates: %d (inside=%d, nearSurface=%d)\n', ...
    size(candMM, 1), sum(candInside), sum(candNear));
fprintf('[AutoElec] Surface support fraction: %.2f\n', surfaceSupport);

if isempty(candMM)
    WC = zeros(0, 3);
    T = (TMax + TMin) / 2;
    fprintf('[AutoElec] No candidates survived peak extraction.\n');
    return;
end

chains = proposeChainsCountGuided(candMM, candConf, candInside, candNear);
fprintf('[AutoElec] Shaft hypotheses: %d\n', numel(chains));

if useCountGuidance
    shafts = selectAndFitShaftsCountGuided(chains, expected.counts, Img, blob, MRInfo.mat);
elseif useChainOptimization
    shafts = selectAndFitShaftsChainOptimization(chains, Img, blob, MRInfo.mat, ...
        ProjSurfRaw, size(candMM, 1), surfaceSupport);
else
    shafts = selectAndFitShaftsAutoCount(chains, Img, blob, MRInfo.mat, ...
        ProjSurfRaw, size(candMM, 1), surfaceSupport);
end
if isempty(shafts)
    WC = zeros(0, 3);
    T = (TMax + TMin) / 2;
    fprintf('[AutoElec] No shafts survived fitting.\n');
    return;
end

finalMM = zeros(0, 3);
for i = 1:numel(shafts)
    finalMM = [finalMM; shafts(i).contactsMM]; %#ok<AGROW>
end

vox = mmToVoxCountGuided(finalMM, MRInfo.mat);
imSz = getOutputImageSizeCountGuided(app, Img);
margin = getDrawSafetyMarginCountGuided(app, imSz);
lo = 1 + margin;
hi = imSz - margin;
hi = max(hi, lo);
keep = vox(:,1) >= lo(1) & vox(:,1) <= hi(1) & ...
       vox(:,2) >= lo(2) & vox(:,2) <= hi(2) & ...
       vox(:,3) >= lo(3) & vox(:,3) <= hi(3);
vox = vox(keep, :);
WC = round(vox);
if ~isempty(WC)
    WC = unique(WC, 'rows', 'stable');
end
T = (TMax + TMin) / 2;

if size(finalMM, 1) ~= size(WC, 1)
    fprintf('[AutoElec] Bounds filter: kept %d / %d contacts within draw-safe MR bounds\n', ...
        size(WC, 1), size(finalMM, 1));
end

fprintf('[AutoElec] FINAL: %d shafts, %d contacts in %.1f s\n', ...
    numel(shafts), size(WC, 1), toc(StartTime));

% % Keep a lightweight threshold curve for backward-compatible diagnostics.
% if ~isempty(ThrR) && ~isempty(ThrHU)
%     try
%         NumObj = zeros(numel(ThrR), 1);
%         for kt = 1:numel(ThrR)
%             CC = bwconncomp(Img > ThrR(kt), 26);
%             if CC.NumObjects == 0
%                 continue;
%             end
%             NumObj(kt) = CC.NumObjects;
%         end
%         [~, idxD] = min(abs(NumObj - size(WC,1)));
%         diagPlotCountGuided(app, ThrHU, NumObj, idxD, toc(StartTime));
%     catch
%     end
% end

end

function expected = getExpectedCountsCountGuided(app)
expected.names = cell(0, 1);
expected.counts = zeros(0, 1);

mapRaw = [];
depthIdx = [];
if isfield(app, 'ElecMapRaw') && ~isempty(app.ElecMapRaw)
    mapRaw = app.ElecMapRaw;
end
if isfield(app, 'DepthElecRaw') && ~isempty(app.DepthElecRaw)
    depthIdx = app.DepthElecRaw(:);
end

if (isempty(mapRaw) || isempty(depthIdx)) && isfield(app, 'ElectFile') && ~isempty(app.ElectFile)
    try
        loaded = load(app.ElectFile);
        if isempty(mapRaw) && isfield(loaded, 'ElecMapRaw')
            mapRaw = loaded.ElecMapRaw;
        end
        if isempty(depthIdx) && isfield(loaded, 'DepthElecRaw')
            depthIdx = loaded.DepthElecRaw(:);
        end
    catch
    end
end

if isempty(mapRaw) || isempty(depthIdx)
    return;
end
if ~iscell(mapRaw) || size(mapRaw, 2) < 1
    return;
end

depthIdx = depthIdx(depthIdx >= 1 & depthIdx <= size(mapRaw, 1));
if isempty(depthIdx)
    return;
end

labels = mapRaw(depthIdx, 1);
labels = labels(:);
labels = labels(~cellfun(@isempty, labels));
if isempty(labels)
    return;
end

names = regexprep(labels, '\d+$', '');
if isempty(names)
    return;
end

[uNames, ~, ic] = unique(names, 'stable');
counts = accumarray(ic, 1);
[counts, ord] = sort(counts(:), 'descend');
expected.names = uNames(ord);
expected.counts = counts(:);
end

function targetCount = candidatePeakTargetCountGuided(imSz)
nVox = prod(double(imSz));
targetCount = round(nVox / 70000.0);
targetCount = min(max(targetCount, 180), 320);
end

function imSz = getOutputImageSizeCountGuided(app, Img)
imSz = size(Img);
try
    if isfield(app, 'ElecImgProjRaw') && ~isempty(app.ElecImgProjRaw)
        imSz = size(app.ElecImgProjRaw);
    elseif isfield(app, 'GrayImg') && ~isempty(app.GrayImg)
        imSz = size(app.GrayImg);
    elseif isfield(app, 'MRImg') && ~isempty(app.MRImg)
        imSz = size(app.MRImg);
    elseif isfield(app, 'MRInfo') && isfield(app.MRInfo, 'dim') && numel(app.MRInfo.dim) >= 3
        imSz = double(app.MRInfo.dim(1:3));
    end
catch
    imSz = size(Img);
end
imSz = double(imSz(:)');
if numel(imSz) < 3
    imSz(3) = 1;
end
end

function margin = getDrawSafetyMarginCountGuided(app, imSz)
elecRadMm = 2.0;
try
    if isprop(app, 'ElecRadEditH') && ~isempty(app.ElecRadEditH) && ...
            isprop(app.ElecRadEditH, 'Value') && ~isempty(app.ElecRadEditH.Value)
        tmp = str2double(app.ElecRadEditH.Value);
        if isfinite(tmp) && tmp > 0
            elecRadMm = double(tmp);
        end
    elseif isfield(app, 'ElecRadEditH') && ~isempty(app.ElecRadEditH) && ...
            isfield(app.ElecRadEditH, 'Value') && ~isempty(app.ElecRadEditH.Value)
        tmp = str2double(app.ElecRadEditH.Value);
        if isfinite(tmp) && tmp > 0
            elecRadMm = double(tmp);
        end
    end
catch
end

scale = [];
try
    if isfield(app, 'XYZScale') && ~isempty(app.XYZScale)
        scale = double(app.XYZScale(:)');
    end
catch
end
if isempty(scale) || numel(scale) < 3 || any(~isfinite(scale)) || any(scale <= 0)
    scale = [0.7 0.7 0.7];
end
margin = ceil(3 * elecRadMm ./ scale);
margin = max(margin, 0);
margin = min(margin, max(floor((double(imSz(:)') - 1) ./ 2), 0));
end

function cap = candidateCapCountGuided(imSz)
nVox = prod(double(imSz));
cap = round(nVox / 18000.0);
cap = min(max(cap, 700), 1400);
end

function blob = computeBlobScoreMapCountGuided(Img, voxDim)
sigmasMm = [0.4 0.7 1.0 1.3];
blob = zeros(size(Img), 'single');
for s = 1:numel(sigmasMm)
    sigmaMm = sigmasMm(s);
    sigmaVox = sigmaMm ./ voxDim(:)';
    g1 = imgaussfilt3(Img, sigmaVox);
    g2 = imgaussfilt3(Img, sigmaVox * 1.6);
    blob = max(blob, -(g2 - g1) * sigmaMm^2);
end
pos = blob(blob > 0);
if ~isempty(pos)
    scl = prctile(pos, 99.5);
    if scl > 0
        blob = blob ./ scl;
    end
end
end

function [candMM, candConf, blobScore, intScore] = ...
    extractLocalMaxCandidateArraysCountGuided(Img, blob, mat)

v = Img(isfinite(Img));
if isempty(v)
    candMM = zeros(0, 3);
    candConf = zeros(0, 1);
    blobScore = zeros(0, 1);
    intScore = zeros(0, 1);
    return;
end

pct = prctile(v, [95 99 99.5 99.8]);
p95 = pct(1);
p99_5 = pct(3);
p99_8 = pct(4);

iScore = clipCountGuided((Img - p95) ./ max(p99_5 - p95, 1e-6), 0, 2.5);
combo = iScore .* max(blob, 0.15);
comboPos = combo(combo > 0);
if isempty(comboPos)
    candMM = zeros(0, 3);
    candConf = zeros(0, 1);
    blobScore = zeros(0, 1);
    intScore = zeros(0, 1);
    return;
end

mx = imdilate(combo, true(3, 3, 3));
peakBase = (combo >= mx) & (Img >= p95);
peakVals = combo(peakBase);
targetMinPeaks = candidatePeakTargetCountGuided(size(Img));
comboThr = max(prctile(comboPos, 99.1), 0.12);
for pctThr = [99.6 99.4 99.2 99.0 98.8 98.6 98.4 98.2 98.0]
    thr = max(prctile(comboPos, pctThr), 0.12);
    comboThr = thr;
    if nnz(peakVals >= thr) >= targetMinPeaks
        break;
    end
end
peakMask = peakBase & (combo >= comboThr);
[ii, jj, kk] = ind2sub(size(Img), find(peakMask));
nP = numel(ii);
if nP == 0
    candMM = zeros(0, 3);
    candConf = zeros(0, 1);
    blobScore = zeros(0, 1);
    intScore = zeros(0, 1);
    return;
end

vox = zeros(nP, 3);
blobScore = zeros(nP, 1);
intScore = zeros(nP, 1);
candConf = zeros(nP, 1);
imSz = size(Img);

for p = 1:nP
    lo = max([ii(p)-1, jj(p)-1, kk(p)-1], [1 1 1]);
    hi = min([ii(p)+1, jj(p)+1, kk(p)+1], imSz);
    [ri, cj, sk] = ndgrid(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3));
    localCombo = combo(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3));
    pts = [ri(:), cj(:), sk(:)];
    w = max(localCombo(:), 0) + 1e-6;
    vox(p, :) = sum(bsxfun(@times, pts, w), 1) ./ sum(w);

    bLoc = sampleVolumeTrilinearCountGuided(blob, vox(p, :));
    iLoc = sampleVolumeTrilinearCountGuided(Img, vox(p, :));
    blobScore(p) = clipCountGuided(bLoc, 0, 2.0);
    intScore(p) = clipCountGuided((iLoc - p95) ./ max(p99_8 - p95, 1e-6), 0, 2.0);
    candConf(p) = 0.55 * intScore(p) + 0.45 * blobScore(p);
end

candMM = voxToMMCountGuided(vox, mat);
[candMM, candConf, blobScore, intScore] = dedupPointsCountGuided(candMM, candConf, blobScore, intScore, 1.6);
end

function ratio = componentAxisRatioCountGuided(pointsMM)
ratio = 1.0;
if size(pointsMM, 1) < 4
    return;
end
mu = mean(pointsMM, 1);
[~, S, ~] = svd(pointsMM - mu, 'econ');
s = diag(S);
if isempty(s)
    return;
end
ratio = s(1) ./ max(s(end), 1e-6);
end

function [candMM, candConf, blobScore, intScore] = mergeThresholdCandidatesCountGuided( ...
    pointsMM, intensityScores, blobScores, strongMask, thresholdCount)

candMM = zeros(0, 3);
candConf = zeros(0, 1);
blobScore = zeros(0, 1);
intScore = zeros(0, 1);
if isempty(pointsMM)
    return;
end

maxRawCandidates = 3000;
clusterRadiusMm = 1.5;
baseScore = intensityScores .* max(blobScores, 0.01);
if size(pointsMM, 1) > maxRawCandidates
    keepMask = logical(strongMask(:));
    nRemaining = maxRawCandidates - nnz(keepMask);
    if nRemaining > 0
        extras = find(~keepMask);
        [~, ord] = sort(baseScore(extras), 'descend');
        keepMask(extras(ord(1:min(nRemaining, numel(ord))))) = true;
    end
    pointsMM = pointsMM(keepMask, :);
    intensityScores = intensityScores(keepMask);
    blobScores = blobScores(keepMask);
end

if size(pointsMM, 1) <= 1
    clusterLabels = ones(size(pointsMM, 1), 1);
else
    try
        Z = linkage(pdist(pointsMM), 'average');
        clusterLabels = cluster(Z, 'cutoff', clusterRadiusMm, 'criterion', 'distance');
    catch
        adj = pdist2(pointsMM, pointsMM) <= clusterRadiusMm;
        G = graph(adj);
        clusterLabels = conncomp(G)';
    end
end

uLabels = unique(clusterLabels(:))';
for clusterId = uLabels
    members = clusterLabels == clusterId;
    pts = pointsMM(members, :);
    eta = intensityScores(members);
    blobVals = blobScores(members);
    weights = eta .* max(blobVals, 0.01) + 1e-6;
    mergedPt = sum(bsxfun(@times, pts, weights), 1) ./ sum(weights);
    mergedInt = max(eta);
    mergedBlob = max(blobVals);
    persistence = nnz(members) ./ max(thresholdCount, 1);

    candMM(end+1, :) = mergedPt; %#ok<AGROW>
    candConf(end+1, 1) = 0.20 * persistence + ...
        0.45 * clipCountGuided(mergedInt, 0, 2.5) + ...
        0.35 * clipCountGuided(mergedBlob, 0, 2.0); %#ok<AGROW>
    blobScore(end+1, 1) = clipCountGuided(mergedBlob, 0, 2.0); %#ok<AGROW>
    intScore(end+1, 1) = clipCountGuided(mergedInt, 0, 2.5); %#ok<AGROW>
end
end

function voxDim = voxelSizesFromMatCountGuided(mat)
A = double(mat(1:3, 1:3));
voxDim = sqrt(sum(A.^2, 1));
voxDim = max(voxDim(:)', 1e-6);
if numel(voxDim) < 3
    voxDim(3) = 1;
end
end

function [candMM, candConf, blobScore, intScore] = ...
    extractComponentCandidateArraysCountGuided(Img, blob, mat, ProjSurfRaw)

valid = Img(isfinite(Img) & Img > 0);
if isempty(valid)
    candMM = zeros(0, 3);
    candConf = zeros(0, 1);
    blobScore = zeros(0, 1);
    intScore = zeros(0, 1);
    return;
end

pct = prctile(valid, [50 99.9]);
p50 = pct(1);
p99_9 = pct(2);
eta = max(0, (Img - p50) ./ max(p99_9 - p50, 1e-6));
combo = eta .* max(blob, 0.1);

voxDim = voxelSizesFromMatCountGuided(mat);
voxelVolumeMm3 = prod(voxDim);
nMin = max(2, floor((0.75 * 0.15) ./ max(voxelVolumeMm3, 1e-6)));
nMax = max(60, ceil((0.75 * 45.0) ./ max(voxelVolumeMm3, 1e-6)));
thresholds = linspace(0.05, 0.90, 18);

rawPoints = zeros(0, 3);
rawIntScore = zeros(0, 1);
rawBlobScore = zeros(0, 1);
rawStrong = false(0, 1);

for t = 1:numel(thresholds)
    CC = bwconncomp(combo >= thresholds(t), 26);
    if CC.NumObjects == 0
        continue;
    end
    for labelId = 1:CC.NumObjects
        linIdx = CC.PixelIdxList{labelId};
        if isempty(linIdx)
            continue;
        end

        [ii, jj, kk] = ind2sub(size(Img), linIdx);
        compVox = [ii(:) jj(:) kk(:)];
        compIntensity = double(Img(linIdx));
        meanIntensity = mean(compIntensity);
        nVox = numel(linIdx);
        tinyBright = nVox >= 1 && nVox <= 3 && meanIntensity >= p99_9;
        if ~((nVox >= nMin && nVox <= nMax) || tinyBright)
            continue;
        end

        compMM = voxToMMCountGuided(compVox, mat);
        if componentAxisRatioCountGuided(compMM) > 8.0
            continue;
        end

        weights = max(compIntensity - p50, 0) + 1e-6;
        centroidVox = sum(bsxfun(@times, compVox, weights), 1) ./ sum(weights);
        centroidMM = voxToMMCountGuided(centroidVox, mat);
        isInside = LeG_intriangulation(ProjSurfRaw.vertices, ProjSurfRaw.faces, centroidMM) == 1;
        if ~isInside && meanIntensity < p99_9
            continue;
        end

        rawPoints(end+1, :) = centroidMM; %#ok<AGROW>
        rawBlobScore(end+1, 1) = sampleVolumeTrilinearCountGuided(blob, centroidVox); %#ok<AGROW>
        rawIntScore(end+1, 1) = sampleVolumeTrilinearCountGuided(eta, centroidVox); %#ok<AGROW>
        rawStrong(end+1, 1) = meanIntensity >= p99_9; %#ok<AGROW>
    end
end

if isempty(rawPoints)
    candMM = zeros(0, 3);
    candConf = zeros(0, 1);
    blobScore = zeros(0, 1);
    intScore = zeros(0, 1);
    return;
end

[candMM, candConf, blobScore, intScore] = mergeThresholdCandidatesCountGuided( ...
    rawPoints, rawIntScore, rawBlobScore, rawStrong, numel(thresholds));
end

function [candMM, candConf, inside, nearSurface, blobScore, intScore] = ...
    extractCandidatesCountGuided(Img, blob, mat, ProjSurfRaw)

v = Img(isfinite(Img));
if isempty(v)
    candMM = zeros(0, 3);
    candConf = zeros(0, 1);
    inside = false(0, 1);
    nearSurface = false(0, 1);
    blobScore = zeros(0, 1);
    intScore = zeros(0, 1);
    return;
end

targetMinPeaks = candidatePeakTargetCountGuided(size(Img));
[candMM, candConf, blobScore, intScore] = ...
    extractComponentCandidateArraysCountGuided(Img, blob, mat, ProjSurfRaw);
[peakMM, peakConf, peakBlob, peakInt] = ...
    extractLocalMaxCandidateArraysCountGuided(Img, blob, mat);

if isempty(candMM)
    candMM = peakMM;
    candConf = peakConf;
    blobScore = peakBlob;
    intScore = peakInt;
else
    rescueThreshold = max(80, round(0.55 * targetMinPeaks));
    if size(candMM, 1) < rescueThreshold && ~isempty(peakMM)
        [candMM, candConf, blobScore, intScore] = dedupPointsCountGuided( ...
            [candMM; peakMM], ...
            [candConf; peakConf], ...
            [blobScore; peakBlob], ...
            [intScore; peakInt], ...
            1.4);
    end
end

if isempty(candMM)
    candMM = peakMM;
    candConf = peakConf;
    blobScore = peakBlob;
    intScore = peakInt;
end

if isempty(candMM)
    inside = false(0, 1);
    nearSurface = false(0, 1);
    return;
end

inside = LeG_intriangulation(ProjSurfRaw.vertices, ProjSurfRaw.faces, candMM) == 1;
dist = nearestSurfaceDistanceCountGuided(candMM, ProjSurfRaw.vertices);
nearSurface = dist <= 6.0;
supportFrac = mean(inside | nearSurface);
if numel(dist) >= 120 && supportFrac < 0.30
    nearSurface = dist <= 10.0;
    supportFrac = mean(inside | nearSurface);
elseif numel(dist) >= 120 && supportFrac < 0.40
    nearSurface = dist <= 8.0;
    supportFrac = mean(inside | nearSurface);
end

hiConfThr = 1.05;
loConfThr = 0.90;
if supportFrac < 0.25
    hiConfThr = 0.95;
    loConfThr = 0.78;
elseif supportFrac < 0.35
    hiConfThr = 1.00;
    loConfThr = 0.84;
end

usable = inside | nearSurface | candConf >= hiConfThr;
minKeep = min(max(160, round(0.75 * targetMinPeaks)), numel(candConf));
if nnz(usable) < minKeep
    usable = inside | nearSurface | candConf >= loConfThr;
end
candMM = candMM(usable, :);
candConf = candConf(usable);
inside = inside(usable);
nearSurface = nearSurface(usable);
blobScore = blobScore(usable);
intScore = intScore(usable);

[~, ord] = sort(candConf, 'descend');
maxCandidates = candidateCapCountGuided(size(Img));
ord = ord(1:min(maxCandidates, numel(ord)));
candMM = candMM(ord, :);
candConf = candConf(ord);
inside = inside(ord);
nearSurface = nearSurface(ord);
blobScore = blobScore(ord);
intScore = intScore(ord);
end

function d = nearestSurfaceDistanceCountGuided(points, vertices)
if isempty(points)
    d = zeros(0, 1);
    return;
end
if size(vertices, 1) > 12000
    step = ceil(size(vertices, 1) / 12000);
    vertices = vertices(1:step:end, :);
end
d = min(pdist2(points, vertices), [], 2);
end

function [pts, conf, blobScore, intScore] = dedupPointsCountGuided(pts, conf, blobScore, intScore, radiusMm)
if size(pts, 1) < 2
    return;
end
[~, ord] = sort(conf, 'descend');
pts = pts(ord, :);
conf = conf(ord);
blobScore = blobScore(ord);
intScore = intScore(ord);

keep = true(size(conf));
for i = 1:size(pts, 1)
    if ~keep(i)
        continue;
    end
    if i == size(pts, 1)
        continue;
    end
    d = sqrt(sum((pts(i+1:end, :) - pts(i, :)).^2, 2));
    keep(i + find(d < radiusMm)) = false;
end

pts = pts(keep, :);
conf = conf(keep);
blobScore = blobScore(keep);
intScore = intScore(keep);
end

function chains = proposeChainsCountGuided(pointsMM, conf, inside, nearSurface)
n = size(pointsMM, 1);
chains = struct('idx', {}, 'pointsMM', {}, 'scores', {}, ...
    'spacing', {}, 'span', {}, 'regularity', {}, 'brainSupport', {}, ...
    'insideFrac', {}, 'nearFrac', {}, 'score', {});
if n < 3
    return;
end

[nbrIdx, nbrDist] = buildNeighborMapCountGuided(pointsMM, 1.8, 6.0);
tripSeeds = zeros(0, 4);
pairSeeds = zeros(0, 3);

for b = 1:n
    nb = nbrIdx{b};
    db = nbrDist{b};
    if isempty(nb)
        continue;
    end
    for k = 1:numel(nb)
        a = nb(k);
        if a < b
            pairSeeds(end+1, :) = [pairSeedScoreCountGuided(conf, a, b, db(k)), a, b]; %#ok<AGROW>
        end
    end
    if numel(nb) < 2
        continue;
    end
    for i = 1:numel(nb)-1
        a = nb(i);
        va = pointsMM(a, :) - pointsMM(b, :);
        da = norm(va);
        if da < 1e-6
            continue;
        end
        for j = i+1:numel(nb)
            c = nb(j);
            vc = pointsMM(c, :) - pointsMM(b, :);
            dc = norm(vc);
            if dc < 1e-6 || abs(da - dc) > 1.5
                continue;
            end
            col = dot(va, vc) / (da * dc);
            if col > -0.90
                continue;
            end
            score = conf(a) + conf(b) + conf(c) + (-col);
            tripSeeds(end+1, :) = [score, a, b, c]; %#ok<AGROW>
        end
    end
end

if ~isempty(tripSeeds)
    [~, ord] = sort(tripSeeds(:, 1), 'descend');
    maxSeeds = min(max(400, round(3.0 * n)), 2400);
    tripSeeds = tripSeeds(ord(1:min(size(tripSeeds,1), maxSeeds)), :);
end
if ~isempty(pairSeeds)
    [~, ord] = sort(pairSeeds(:, 1), 'descend');
    maxPairSeeds = min(max(600, round(4.0 * n)), 3200);
    pairSeeds = pairSeeds(ord(1:min(size(pairSeeds,1), maxPairSeeds)), :);
end

keys = cell(0, 1);
for s = 1:size(tripSeeds, 1)
    chain = growChainFromSeedCountGuided(tripSeeds(s, 2:4), pointsMM, conf, inside, nearSurface);
    if isempty(chain)
        continue;
    end
    key = sprintf('%d_', sort(chain.idx));
    if any(strcmp(keys, key))
        continue;
    end
    keys{end+1,1} = key; %#ok<AGROW>
    chains(end+1) = chain; %#ok<AGROW>
end
for s = 1:size(pairSeeds, 1)
    chain = growChainFromPairCountGuided(pairSeeds(s, 2:3), pointsMM, conf, inside, nearSurface);
    if isempty(chain)
        continue;
    end
    key = sprintf('%d_', sort(chain.idx));
    if any(strcmp(keys, key))
        continue;
    end
    keys{end+1,1} = key; %#ok<AGROW>
    chains(end+1) = chain; %#ok<AGROW>
end

if ~isempty(chains)
    [~, ord] = sort([chains.score], 'descend');
    chains = chains(ord);
end
end

function [nbrIdx, nbrDist] = buildNeighborMapCountGuided(pointsMM, minStep, maxStep)
n = size(pointsMM, 1);
nbrIdx = cell(n, 1);
nbrDist = cell(n, 1);
if n < 2
    return;
end

D = pdist2(pointsMM, pointsMM);
for i = 1:n
    mask = D(i, :) >= minStep & D(i, :) <= maxStep;
    mask(i) = false;
    nbrIdx{i} = find(mask);
    nbrDist{i} = D(i, mask);
end
end

function score = pairSeedScoreCountGuided(conf, a, b, d)
score = conf(a) + conf(b) - 0.25 * abs(d - 3.5);
end

function chain = growChainFromSeedCountGuided(seedIdx, pointsMM, conf, inside, nearSurface)
chain = [];
seedPts = pointsMM(seedIdx, :);
mu = mean(seedPts, 1);
[~, ~, V] = svd(seedPts - mu, 'econ');
axis1 = V(:, 1);
t = (seedPts - mu) * axis1;
[~, ord] = sort(t);
idx = seedIdx(ord);
spacing = estimateSpacingCountGuided(pointsMM(idx, :));
if spacing < 2.0 || spacing > 5.5
    spacing = 3.5;
end

idx = extendChainCountGuided(idx, pointsMM, conf, spacing, true);
idx = extendChainCountGuided(idx, pointsMM, conf, spacing, false);
idx = removeDuplicateIndicesCountGuided(idx);
if numel(idx) < 3
    return;
end

pts = pointsMM(idx, :);
scores = conf(idx);
spacing = estimateSpacingCountGuided(pts);
regularity = spacingRegularityCountGuided(pts, spacing);
span = polylineLengthCountGuided(pts);
insideFrac = mean(inside(idx));
nearFrac = mean(nearSurface(idx));
brainSupport = mean(inside(idx) | nearSurface(idx));
if regularity < 0.22 || span < 6.0
    return;
end

chain.idx = idx(:)';
chain.pointsMM = pts;
chain.scores = scores;
chain.spacing = min(max(spacing, 2.4), 5.0);
chain.span = span;
chain.regularity = regularity;
chain.brainSupport = brainSupport;
chain.insideFrac = insideFrac;
chain.nearFrac = nearFrac;
chain.score = 1.7 * regularity + 0.8 * mean(scores) + 0.24 * numel(idx) + ...
    0.06 * span + 0.45 * brainSupport + 0.12 * insideFrac;
end

function chain = growChainFromPairCountGuided(pairIdx, pointsMM, conf, inside, nearSurface)
chain = [];
a = pairIdx(1);
b = pairIdx(2);
d = norm(pointsMM(b, :) - pointsMM(a, :));
if d < 1.8 || d > 6.2
    return;
end

idx = extendChainCountGuided([a b], pointsMM, conf, d, true);
idx = extendChainCountGuided(idx, pointsMM, conf, d, false);
idx = removeDuplicateIndicesCountGuided(idx);
if numel(idx) < 3
    return;
end

pts = pointsMM(idx, :);
scores = conf(idx);
spacing = estimateSpacingCountGuided(pts);
regularity = spacingRegularityCountGuided(pts, spacing);
span = polylineLengthCountGuided(pts);
insideFrac = mean(inside(idx));
nearFrac = mean(nearSurface(idx));
brainSupport = mean(inside(idx) | nearSurface(idx));
if regularity < 0.18 || span < 6.0
    return;
end

chain.idx = idx(:)';
chain.pointsMM = pts;
chain.scores = scores;
chain.spacing = min(max(spacing, 2.4), 5.0);
chain.span = span;
chain.regularity = regularity;
chain.brainSupport = brainSupport;
chain.insideFrac = insideFrac;
chain.nearFrac = nearFrac;
chain.score = 1.3 * regularity + 0.8 * mean(scores) + 0.22 * numel(idx) + ...
    0.05 * span + 0.40 * brainSupport + 0.10 * insideFrac;
end

function idx = extendChainCountGuided(idx, pointsMM, conf, spacing, isForward)
maxLen = 22;
idx = idx(:)';
while numel(idx) < maxLen
    if isForward
        work = idx;
    else
        work = fliplr(idx);
    end
    pts = pointsMM(work, :);
    if size(pts, 1) >= 3
        tangent = pts(end, :) - pts(end-2, :);
    else
        tangent = pts(end, :) - pts(end-1, :);
    end
    nrm = norm(tangent);
    if nrm < 1e-6
        break;
    end
    tangent = tangent ./ nrm;
    tip = pts(end, :);

    bestIdx = 0;
    bestScore = -Inf;
    for c = 1:size(pointsMM, 1)
        if any(idx == c)
            continue;
        end
        stepVec = pointsMM(c, :) - tip;
        step = norm(stepVec);
        if step < 1.6 || step > 6.2
            continue;
        end
        align = dot(stepVec ./ step, tangent);
        if align < 0.78
            continue;
        end
        pred = tip + tangent .* spacing;
        predErr = norm(pointsMM(c, :) - pred);
        score = 1.8 * conf(c) + 1.2 * align - 0.55 * predErr - 0.18 * abs(step - spacing);
        if score > bestScore
            bestScore = score;
            bestIdx = c;
        end
    end

    if bestIdx == 0 || bestScore < 0.45
        break;
    end

    if isForward
        idx(end+1) = bestIdx; %#ok<AGROW>
    else
        idx = [bestIdx idx]; %#ok<AGROW>
    end
    spacing = estimateSpacingCountGuided(pointsMM(idx, :));
end
end

function idx = removeDuplicateIndicesCountGuided(idx)
out = zeros(1, 0);
for i = 1:numel(idx)
    if ~any(out == idx(i))
        out(end+1) = idx(i); %#ok<AGROW>
    end
end
idx = out;
end

function spacing = estimateSpacingCountGuided(pointsMM)
spacing = 3.5;
if size(pointsMM, 1) < 2
    return;
end
t = projectAxisCountGuided(pointsMM);
gaps = diff(sort(t));
gaps = gaps(gaps > 1.5 & gaps < 6.5);
if isempty(gaps)
    return;
end
spacing = median(gaps);
end

function reg = spacingRegularityCountGuided(pointsMM, spacing)
reg = 0;
if size(pointsMM, 1) < 3
    return;
end
t = sort(projectAxisCountGuided(pointsMM));
gaps = diff(t);
gaps = gaps(gaps > 1.2 & gaps < 7.0);
if isempty(gaps)
    return;
end
err = abs(gaps - spacing);
reg = mean(exp(-((err ./ 1.0).^2)));
end

function t = projectAxisCountGuided(pointsMM)
mu = mean(pointsMM, 1);
[~, ~, V] = svd(pointsMM - mu, 'econ');
axis1 = V(:, 1);
t = (pointsMM - mu) * axis1;
end

function L = polylineLengthCountGuided(pointsMM)
if size(pointsMM, 1) < 2
    L = 0;
    return;
end
L = sum(sqrt(sum(diff(pointsMM, 1, 1).^2, 2)));
end

function shafts = selectAndFitShaftsAutoCount(chains, Img, blob, mat, ProjSurfRaw, nCandidates, surfaceSupport)
shafts = emptyAutoCountShaftStruct();
if isempty(chains)
    return;
end

nChains = min(numel(chains), 140);
chains = chains(1:nChains);
fits = emptyAutoCountShaftStruct();
for i = 1:nChains
    shaft = fitContactsToChainAutoCount(chains(i), Img, blob, mat, ProjSurfRaw, surfaceSupport);
    if isempty(shaft)
        continue;
    end
    fits(end+1) = shaft; %#ok<AGROW>
end
if isempty(fits)
    return;
end

scores = [fits.score];
scoreThr = max(2.20, prctile(scores, 40) - 0.35);
if ~isfinite(scoreThr)
    scoreThr = 2.20;
end
if surfaceSupport < 0.35
    scoreThr = scoreThr - 0.20;
elseif surfaceSupport < 0.45
    scoreThr = scoreThr - 0.10;
elseif surfaceSupport < 0.55
    scoreThr = scoreThr - 0.05;
end

minSupport = 0.32;
minEvidence = 0.32;
minBrainSupport = 0.30;
maxOutsideFrac = 0.45;
brainEvidenceBypass = 0.75;
if surfaceSupport < 0.25
    minBrainSupport = 0.06;
    maxOutsideFrac = 0.78;
    brainEvidenceBypass = 0.55;
elseif surfaceSupport < 0.35
    minBrainSupport = 0.12;
    maxOutsideFrac = 0.68;
    brainEvidenceBypass = 0.60;
elseif surfaceSupport < 0.45
    minBrainSupport = 0.20;
    maxOutsideFrac = 0.56;
    brainEvidenceBypass = 0.68;
end

[~, ord] = sort([fits.score], 'descend');
fits = fits(ord);

maxIdx = 0;
for i = 1:numel(chains)
    if ~isempty(chains(i).idx)
        maxIdx = max(maxIdx, max(chains(i).idx));
    end
end
usedMask = false(1, maxIdx);
selected = emptyAutoCountShaftStruct();
totalContacts = 0;
maxTotalContacts = max(320, 2 * max(1, round(double(nCandidates))));

for i = 1:numel(fits)
    shaft = fits(i);
    overlap = 0;
    if ~isempty(shaft.seedIdx) && ~isempty(usedMask)
        overlap = mean(usedMask(shaft.seedIdx));
    end
    if overlap > 0.58
        continue;
    end
    if shaft.score < scoreThr
        continue;
    end
    if shaft.supportScore < minSupport || shaft.evidenceScore < minEvidence
        continue;
    end
    if shaft.brainSupport < minBrainSupport && shaft.evidenceScore < brainEvidenceBypass
        continue;
    end
    if shaft.outsideFrac > maxOutsideFrac && shaft.evidenceScore < 0.95
        continue;
    end
    if isDuplicateShaftAutoCountGuided(selected, shaft)
        continue;
    end
    if totalContacts + size(shaft.contactsMM, 1) > maxTotalContacts
        continue;
    end
    selected(end+1) = shaft; %#ok<AGROW>
    if ~isempty(shaft.seedIdx)
        usedMask(shaft.seedIdx) = true;
    end
    totalContacts = totalContacts + size(shaft.contactsMM, 1);
end

selected = addRelaxedExtraShaftsAutoCount(selected, fits, scoreThr, maxTotalContacts, ...
    minSupport, minEvidence, minBrainSupport, maxOutsideFrac, ...
    brainEvidenceBypass, surfaceSupport);

if isempty(selected)
    shafts = selected;
    return;
end

if surfaceSupport < 0.55
    selected = mergeCollinearShaftsAutoCount(selected, false);
    selected = consolidateShaftGroupsAutoCount(selected, Img, blob, mat, ...
        ProjSurfRaw, surfaceSupport, false);
end
selected = inflateSelectedShaftsAutoCount(selected, Img, blob, mat, ProjSurfRaw, surfaceSupport);
strictPostMerge = surfaceSupport >= 0.60;
selected = mergeCollinearShaftsAutoCount(selected, strictPostMerge);
selected = consolidateShaftGroupsAutoCount(selected, Img, blob, mat, ...
    ProjSurfRaw, surfaceSupport, strictPostMerge);
nContactsBeforeTrim = 0;
for i = 1:numel(selected)
    nContactsBeforeTrim = nContactsBeforeTrim + size(selected(i).contactsMM, 1);
end
selected = trimOutsideShaftEndsAutoCount(selected, Img, blob, mat, ProjSurfRaw, surfaceSupport);
nContactsAfterTrim = 0;
for i = 1:numel(selected)
    nContactsAfterTrim = nContactsAfterTrim + size(selected(i).contactsMM, 1);
end
if nContactsAfterTrim ~= nContactsBeforeTrim
    fprintf('[AutoElec] Outside trim: %d -> %d contacts after removing extracranial shaft ends\n', ...
        nContactsBeforeTrim, nContactsAfterTrim);
end
selected = dropRedundantFragmentsAutoCount(selected);
selected = rescueIntracranialFitsAutoCount(selected, fits, scoreThr, surfaceSupport);
selected = mergeCollinearShaftsAutoCount(selected, strictPostMerge);
selected = consolidateShaftGroupsAutoCount(selected, Img, blob, mat, ...
    ProjSurfRaw, surfaceSupport, strictPostMerge);
selected = trimOutsideShaftEndsAutoCount(selected, Img, blob, mat, ProjSurfRaw, surfaceSupport);
selected = dropRedundantFragmentsAutoCount(selected);

pruned = emptyAutoCountShaftStruct();
for i = 1:numel(selected)
    shaft = selected(i);
    spanMm = polylineLengthCountGuided(shaft.contactsMM);
    if shaft.expectedContacts >= 5 && shaft.insideFrac <= 0.02 && ...
            shaft.outsideFrac >= 0.65 && shaft.evidenceScore < 0.92
        continue;
    end
    if shaft.score < 5.6 && shaft.supportScore < 0.55 && shaft.evidenceScore < 0.90
        continue;
    end
    if shaft.expectedContacts <= 2 && (spanMm < 8.0 || shaft.evidenceScore < 0.90)
        continue;
    end
    if shaft.expectedContacts <= 3 && spanMm < 10.0 && shaft.score < scoreThr + 0.35
        continue;
    end
    if shaft.expectedContacts <= 4 && spanMm < 12.0 && shaft.brainSupport < 0.18 && ...
            shaft.evidenceScore < 0.75
        continue;
    end
    pruned(end+1) = shaft; %#ok<AGROW>
end

if isempty(pruned)
    shafts = pruned;
    return;
end
[~, ord] = sort([pruned.score], 'descend');
shafts = pruned(ord);
end

function shafts = selectAndFitShaftsChainOptimization(chains, Img, blob, mat, ProjSurfRaw, nCandidates, surfaceSupport)
shafts = emptyAutoCountShaftStruct();
if isempty(chains)
    return;
end

nChains = min(numel(chains), 180);
chains = chains(1:nChains);
fits = emptyAutoCountShaftStruct();
for i = 1:nChains
    shaft = fitContactsToChainAutoCount(chains(i), Img, blob, mat, ProjSurfRaw, surfaceSupport);
    if isempty(shaft)
        continue;
    end
    fits(end+1) = shaft; %#ok<AGROW>
end
if isempty(fits)
    return;
end

scores = [fits.score];
scoreThr = max(2.20, prctile(scores, 40) - 0.35);
if ~isfinite(scoreThr)
    scoreThr = 2.20;
end
if surfaceSupport < 0.35
    scoreThr = scoreThr - 0.20;
elseif surfaceSupport < 0.45
    scoreThr = scoreThr - 0.10;
elseif surfaceSupport < 0.55
    scoreThr = scoreThr - 0.05;
end

minSupport = 0.32;
minEvidence = 0.32;
minBrainSupport = 0.30;
maxOutsideFrac = 0.45;
brainEvidenceBypass = 0.75;
if surfaceSupport < 0.25
    minBrainSupport = 0.06;
    maxOutsideFrac = 0.78;
    brainEvidenceBypass = 0.55;
elseif surfaceSupport < 0.35
    minBrainSupport = 0.12;
    maxOutsideFrac = 0.68;
    brainEvidenceBypass = 0.60;
elseif surfaceSupport < 0.45
    minBrainSupport = 0.20;
    maxOutsideFrac = 0.56;
    brainEvidenceBypass = 0.68;
end

[~, ord] = sort([fits.score], 'descend');
fits = fits(ord);
candidatePool = emptyAutoCountShaftStruct();
for i = 1:numel(fits)
    shaft = fits(i);
    if shaft.score < scoreThr - 0.45 && shaft.evidenceScore < 0.82
        continue;
    end
    if shaft.supportScore < minSupport || shaft.evidenceScore < minEvidence
        continue;
    end
    if shaft.brainSupport < minBrainSupport && shaft.evidenceScore < brainEvidenceBypass
        continue;
    end
    if shaft.outsideFrac > maxOutsideFrac + 0.12 && shaft.evidenceScore < 0.95
        continue;
    end
    candidatePool(end+1) = shaft; %#ok<AGROW>
end
if isempty(candidatePool)
    return;
end

maxTotalContacts = max(320, 2 * max(1, round(double(nCandidates))));
selected = optimizeShaftSubsetAutoCount(candidatePool, maxTotalContacts, surfaceSupport, true, false, 256);
selected = addRelaxedExtraShaftsAutoCount(selected, fits, scoreThr, maxTotalContacts, ...
    minSupport, minEvidence, minBrainSupport, maxOutsideFrac, ...
    brainEvidenceBypass, surfaceSupport);
if isempty(selected)
    return;
end

if surfaceSupport < 0.55
    selected = mergeCollinearShaftsAutoCount(selected, false);
    selected = consolidateShaftGroupsAutoCount(selected, Img, blob, mat, ...
        ProjSurfRaw, surfaceSupport, false);
end
selected = inflateSelectedShaftsAutoCount(selected, Img, blob, mat, ProjSurfRaw, surfaceSupport);
strictPostMerge = surfaceSupport >= 0.60;
selected = mergeCollinearShaftsAutoCount(selected, strictPostMerge);
selected = consolidateShaftGroupsAutoCount(selected, Img, blob, mat, ...
    ProjSurfRaw, surfaceSupport, strictPostMerge);
selected = dropRedundantFragmentsAutoCount(selected);

totalContacts = 0;
for i = 1:numel(selected)
    totalContacts = totalContacts + size(selected(i).contactsMM, 1);
end
meanContacts = 0;
if ~isempty(selected)
    meanContacts = totalContacts ./ numel(selected);
end
if numel(selected) >= 15 && meanContacts >= 10.5
    selected = optimizeShaftSubsetAutoCount(selected, maxTotalContacts, surfaceSupport, false, true, 128);
end

pruned = emptyAutoCountShaftStruct();
for i = 1:numel(selected)
    shaft = selected(i);
    spanMm = polylineLengthCountGuided(shaft.contactsMM);
    if shaft.score < 5.6 && shaft.supportScore < 0.55 && shaft.evidenceScore < 0.90
        continue;
    end
    if shaft.expectedContacts <= 2 && (spanMm < 8.0 || shaft.evidenceScore < 0.90)
        continue;
    end
    if shaft.expectedContacts <= 3 && spanMm < 10.0 && shaft.score < scoreThr + 0.35
        continue;
    end
    if shaft.expectedContacts <= 4 && spanMm < 12.0 && shaft.brainSupport < 0.18 && shaft.evidenceScore < 0.75
        continue;
    end
    pruned(end+1) = shaft; %#ok<AGROW>
end

if isempty(pruned)
    shafts = pruned;
    return;
end
[~, ord] = sort([pruned.score], 'descend');
shafts = pruned(ord);
end

function selected = addRelaxedExtraShaftsAutoCount(selected, fits, scoreThr, maxTotalContacts, ...
    minSupport, minEvidence, minBrainSupport, maxOutsideFrac, brainEvidenceBypass, surfaceSupport)
if isempty(fits)
    return;
end

relaxedScoreThr = scoreThr;
extraSlots = 0;
if surfaceSupport >= 0.85
    relaxedScoreThr = max(4.85, scoreThr - 0.90);
    extraSlots = 5;
elseif surfaceSupport >= 0.70
    relaxedScoreThr = max(5.00, scoreThr - 0.60);
    extraSlots = 3;
elseif surfaceSupport >= 0.55
    relaxedScoreThr = max(5.10, scoreThr - 0.45);
    extraSlots = 2;
end
if extraSlots <= 0
    return;
end

maxIdx = 0;
for i = 1:numel(fits)
    if ~isempty(fits(i).seedIdx)
        maxIdx = max(maxIdx, max(fits(i).seedIdx));
    end
end
for i = 1:numel(selected)
    if ~isempty(selected(i).seedIdx)
        maxIdx = max(maxIdx, max(selected(i).seedIdx));
    end
end
usedMask = false(1, maxIdx);
totalContacts = 0;
for i = 1:numel(selected)
    if ~isempty(selected(i).seedIdx)
        seedIdx = selected(i).seedIdx;
        seedIdx = seedIdx(seedIdx >= 1 & seedIdx <= numel(usedMask));
        usedMask(seedIdx) = true;
    end
    totalContacts = totalContacts + size(selected(i).contactsMM, 1);
end

for i = 1:numel(fits)
    shaft = fits(i);
    overlap = 0;
    seedIdx = shaft.seedIdx;
    if ~isempty(seedIdx) && ~isempty(usedMask)
        seedIdx = seedIdx(seedIdx >= 1 & seedIdx <= numel(usedMask));
        if ~isempty(seedIdx)
            overlap = mean(usedMask(seedIdx));
        end
    end
    if overlap > 0.58
        continue;
    end
    if shaft.score < relaxedScoreThr
        continue;
    end
    if shaft.supportScore < minSupport || shaft.evidenceScore < minEvidence
        continue;
    end
    if shaft.brainSupport < minBrainSupport && shaft.evidenceScore < brainEvidenceBypass
        continue;
    end
    if shaft.outsideFrac > maxOutsideFrac && shaft.evidenceScore < 0.95
        continue;
    end
    if shaft.expectedContacts < 6 && shaft.evidenceScore < 0.95
        continue;
    end
    if isDuplicateShaftAutoCountGuided(selected, shaft)
        continue;
    end
    if totalContacts + size(shaft.contactsMM, 1) > maxTotalContacts
        continue;
    end

    selected(end+1) = shaft; %#ok<AGROW>
    if ~isempty(seedIdx)
        usedMask(seedIdx) = true;
    end
    totalContacts = totalContacts + size(shaft.contactsMM, 1);
    extraSlots = extraSlots - 1;
    if extraSlots <= 0
        break;
    end
end
end

function shafts = optimizeShaftSubsetAutoCount(shaftsIn, maxTotalContacts, surfaceSupport, useSeedConflicts, allowFragmentConflicts, beamWidth)
shafts = emptyAutoCountShaftStruct();
if isempty(shaftsIn)
    return;
end

util = zeros(numel(shaftsIn), 1);
for i = 1:numel(shaftsIn)
    util(i) = shaftSelectionUtilityAutoCount(shaftsIn(i), surfaceSupport);
end
[~, ord] = sort(util, 'descend');
ord = ord(1:min(numel(ord), 60));
shaftsIn = shaftsIn(ord);
util = util(ord);
nS = numel(shaftsIn);
if nS == 0
    return;
end

confMasks = buildShaftConflictMasksAutoCount(shaftsIn, useSeedConflicts, allowFragmentConflicts);
suffixPos = zeros(nS + 1, 1);
for i = nS:-1:1
    suffixPos(i) = suffixPos(i + 1) + max(0, util(i));
end

stateVal = 0;
stateSel = uint64(0);
stateForb = uint64(0);
stateTot = 0;

for i = 1:nS
    bit = bitshift(uint64(1), i - 1);
    maxStates = numel(stateVal) * 2;
    newVal = zeros(maxStates, 1);
    newSel = zeros(maxStates, 1, 'uint64');
    newForb = zeros(maxStates, 1, 'uint64');
    newTot = zeros(maxStates, 1);
    idx = 0;
    for s = 1:numel(stateVal)
        idx = idx + 1;
        newVal(idx) = stateVal(s);
        newSel(idx) = stateSel(s);
        newForb(idx) = stateForb(s);
        newTot(idx) = stateTot(s);

        if bitand(stateForb(s), bit) ~= 0
            continue;
        end
        contactAdd = shaftsIn(i).expectedContacts;
        if stateTot(s) + contactAdd > maxTotalContacts
            continue;
        end
        idx = idx + 1;
        newVal(idx) = stateVal(s) + util(i);
        newSel(idx) = bitor(stateSel(s), bit);
        newForb(idx) = bitor(stateForb(s), bitor(confMasks(i), bit));
        newTot(idx) = stateTot(s) + contactAdd;
    end

    newVal = newVal(1:idx);
    newSel = newSel(1:idx);
    newForb = newForb(1:idx);
    newTot = newTot(1:idx);
    optimistic = newVal + suffixPos(i + 1);
    sortKey = [optimistic newVal];
    [~, ord] = sortrows(sortKey, [-1 -2]);
    ord = ord(1:min(numel(ord), beamWidth));
    stateVal = newVal(ord);
    stateSel = newSel(ord);
    stateForb = newForb(ord);
    stateTot = newTot(ord);
end

[~, bestIdx] = max(stateVal);
bestMask = stateSel(bestIdx);
for i = 1:nS
    bit = bitshift(uint64(1), i - 1);
    if bitand(bestMask, bit) == 0
        continue;
    end
    shafts(end+1) = shaftsIn(i); %#ok<AGROW>
end
if ~isempty(shafts)
    [~, ord] = sort([shafts.score], 'descend');
    shafts = shafts(ord);
end
end

function util = shaftSelectionUtilityAutoCount(shaft, surfaceSupport)
surfWeight = 0.25;
if surfaceSupport >= 0.45
    surfWeight = 0.45;
end
util = shaft.score + 0.35 * shaft.supportScore + 0.55 * shaft.evidenceScore + ...
    0.10 * min(shaft.expectedContacts, 12) + ...
    surfWeight * (0.90 * shaft.brainSupport + 0.60 * shaft.insideFrac - 1.10 * shaft.outsideFrac);
end

function confMasks = buildShaftConflictMasksAutoCount(shafts, useSeedConflicts, allowFragmentConflicts)
nS = numel(shafts);
confMasks = zeros(nS, 1, 'uint64');
for i = 1:nS-1
    for j = i+1:nS
        if ~shaftsConflictAutoCount(shafts(i), shafts(j), useSeedConflicts, allowFragmentConflicts)
            continue;
        end
        confMasks(i) = bitor(confMasks(i), bitshift(uint64(1), j - 1));
        confMasks(j) = bitor(confMasks(j), bitshift(uint64(1), i - 1));
    end
end
end

function tf = shaftsConflictAutoCount(a, b, useSeedConflicts, allowFragmentConflicts)
tf = false;
if useSeedConflicts && ~isempty(a.seedIdx) && ~isempty(b.seedIdx)
    overlap = numel(intersect(a.seedIdx(:)', b.seedIdx(:)')) ./ ...
        max(1, min(numel(a.seedIdx), numel(b.seedIdx)));
    if overlap > 0.55
        tf = true;
        return;
    end
end
if shaftsDuplicatePairAutoCount(a, b)
    tf = true;
    return;
end
if allowFragmentConflicts && (isRedundantFragmentAutoCount(a, b) || isRedundantFragmentAutoCount(b, a))
    tf = true;
end
end

function tf = shaftsDuplicatePairAutoCount(a, b)
tf = false;
if isempty(a.contactsMM) || isempty(b.contactsMM)
    return;
end
D = pdist2(a.contactsMM, b.contactsMM);
fracA = mean(min(D, [], 2) < 1.6);
fracB = mean(min(D, [], 1) < 1.6);
tf = fracA > 0.60 || fracB > 0.60;
end

function shaft = fitContactsToChainAutoCount(chain, Img, blob, mat, ProjSurfRaw, surfaceSupport)
shaft = [];
pts = chain.pointsMM;
obsScores = chain.scores;
obsS = arcLengthsCountGuided(pts);
spacingGuess = min(max(chain.spacing, 2.6), 4.6);
span = max(obsS(end), spacingGuess);
obsCount = numel(chain.idx);
baseCount = max(3, round(span ./ max(spacingGuess, 1e-3)) + 1);
gapCount = 0;
if numel(obsS) > 1
    gapCount = sum(diff(obsS) > 1.45 * spacingGuess);
end
countMin = max(3, min([obsCount, baseCount]) - 1);
countMax = min(22, max([obsCount + gapCount + 3, baseCount + 5, ceil((span + 12.0) / 2.8) + 1]));
countGrid = countMin:countMax;
if isempty(countGrid)
    countGrid = max(3, obsCount):(max(3, obsCount) + 2);
end

v = Img(isfinite(Img));
pct = prctile(v, [95 99.8]);
p95 = pct(1);
p99_8 = pct(2);

bestScore = -Inf;
bestSupport = 0;
bestEvidence = 0;
bestSpacing = spacingGuess;
bestPts = zeros(0, 3);
bestGridS = zeros(0, 1);
bestCount = 0;

spacingGrid = max(2.5, spacingGuess - 0.6):0.15:min(4.9, spacingGuess + 0.6);
for spacing = spacingGrid
    for nContacts = countGrid
        nominalLen = spacing * max(nContacts - 1, 1);
        startMin = min(obsS(1), span - nominalLen) - min(3.0, spacing);
        startMax = max(obsS(end) - nominalLen, 0.0) + min(4.5, 1.5 * spacing);
        for s0 = startMin:0.5:startMax
            gridS = s0 + spacing * (0:nContacts-1)';
            gridPts = samplePolylineExtrapCountGuided(pts, gridS);
            vox = mmToVoxCountGuided(gridPts, mat);
            ctVals = sampleVolumeTrilinearCountGuided(Img, vox);
            blobVals = sampleVolumeTrilinearCountGuided(blob, vox);
            ctScore = clipCountGuided((ctVals - p95) ./ max(p99_8 - p95, 1e-6), 0, 2.0);
            blobScore = clipCountGuided(blobVals, 0, 2.0);
            evidenceVec = 0.60 * ctScore + 0.40 * blobScore;

            dObs = abs(bsxfun(@minus, gridS(:), obsS(:)'));
            supportObs = mean(exp(-((min(dObs, [], 1) ./ 1.10).^2)));
            supportGrid = mean(exp(-((min(dObs, [], 2) ./ 1.35).^2)));
            support = 0.55 * supportObs + 0.45 * supportGrid;
            evidence = mean(evidenceVec);
            peakFrac = mean(evidenceVec >= 0.50);
            extrapMm = max(0, -min(gridS)) + max(0, max(gridS) - obsS(end));
            extrapPenalty = extrapMm ./ max(nominalLen, spacing);
            inflatePenalty = max(0, nContacts - baseCount - 2);
            overspanPenalty = max(0, nContacts - max(baseCount, obsCount + 2));

            score = 1.70 * supportObs + 0.40 * supportGrid + ...
                1.80 * evidence + 0.90 * peakFrac + 1.00 * chain.regularity + ...
                0.22 * chain.brainSupport + 0.12 * mean(obsScores) + ...
                0.02 * min(nContacts, 10) - 0.08 * abs(spacing - spacingGuess) - ...
                0.28 * extrapPenalty - 0.11 * inflatePenalty - 0.08 * overspanPenalty;
            if score > bestScore
                bestScore = score;
                bestSupport = support;
                bestEvidence = evidence;
                bestSpacing = spacing;
                bestPts = gridPts;
                bestGridS = gridS;
                bestCount = nContacts;
            end
        end
    end
end

if isempty(bestPts)
    return;
end

[supportDist, outsideDist, surfaceWeight, baseEvThr, strongEvThr] = ...
    surfaceDistanceParamsAutoCount(surfaceSupport);

score = bestScore;
support = bestSupport;
evidence = bestEvidence;
spacingMm = bestSpacing;
nContacts = bestCount;
gridS = bestGridS;
gridPts = bestPts;

if surfaceSupport >= 0.65
    maxExtraTotal = 1;
else
    maxExtraTotal = 2;
end
if evidence >= 0.95 && support >= 0.55
    maxExtraTotal = maxExtraTotal + 1;
end
if evidence >= 1.15 && support >= 0.70
    maxExtraTotal = maxExtraTotal + 2;
end
if chain.regularity >= 0.45 && obsCount >= 5
    maxExtraTotal = maxExtraTotal + 1;
end
if obsCount >= 8
    maxExtraTotal = maxExtraTotal + 1;
end
maxExtraCap = 6;
if surfaceSupport >= 0.50
    maxExtraCap = 2;
elseif surfaceSupport >= 0.45
    maxExtraCap = 2;
elseif surfaceSupport >= 0.40
    maxExtraCap = 3;
end
maxExtraTotal = min(max(maxExtraTotal, 1), maxExtraCap);
evRelax = 0;
strongRelax = 0;
if surfaceSupport < 0.40 && evidence >= 1.15
    evRelax = 0.03;
end
if surfaceSupport < 0.35 && evidence >= 1.25
    strongRelax = 0.04;
end
[gridS, gridPts] = extendGridByEvidenceAutoCount(pts, gridS, spacingMm, Img, blob, ...
    mat, p95, p99_8, ProjSurfRaw, supportDist, outsideDist, maxExtraTotal, ...
    max(0.52, baseEvThr - evRelax), max(baseEvThr + 0.18, strongEvThr - strongRelax));

nContacts = numel(gridS);
evidenceVals = contactEvidenceValuesAutoCount(Img, blob, mat, p95, p99_8, gridPts);
evidence = mean(evidenceVals);
peakFrac = mean(evidenceVals >= 0.50);
score = score + 0.10 * min(max(nContacts - bestCount, 0), 6) + 0.20 * peakFrac;

inside = LeG_intriangulation(ProjSurfRaw.vertices, ProjSurfRaw.faces, gridPts) == 1;
dist = nearestSurfaceDistanceCountGuided(gridPts, ProjSurfRaw.vertices);
brainMask = inside | (dist <= supportDist);
brainSupport = mean(brainMask);
outsideFrac = mean(~inside & dist > outsideDist);
insideFrac = mean(inside);
score = score + surfaceWeight * (0.35 * brainSupport - 0.55 * outsideFrac);

if support < 0.30 || evidence < 0.30
    return;
end
if surfaceSupport >= 0.35 && outsideFrac > 0.60 && evidence < 0.80
    return;
end

shaft.expectedContacts = nContacts;
shaft.observedContacts = numel(chain.idx);
shaft.score = score;
shaft.supportScore = support;
shaft.spacingMm = spacingMm;
shaft.contactsMM = gridPts;
shaft.seedIdx = chain.idx;
shaft.brainSupport = brainSupport;
shaft.insideFrac = insideFrac;
shaft.evidenceScore = evidence;
shaft.outsideFrac = outsideFrac;
end

function tf = isDuplicateShaftAutoCountGuided(shafts, shaft)
tf = false;
if isempty(shafts) || isempty(shaft.contactsMM)
    return;
end
for i = 1:numel(shafts)
    other = shafts(i).contactsMM;
    if isempty(other)
        continue;
    end
    D = pdist2(shaft.contactsMM, other);
    fracA = mean(min(D, [], 2) < 1.6);
    fracB = mean(min(D, [], 1) < 1.6);
    if fracA > 0.60 || fracB > 0.60
        tf = true;
        return;
    end
end
end

function vals = contactEvidenceValuesAutoCount(Img, blob, mat, p95, p99_8, pointsMM)
vox = mmToVoxCountGuided(pointsMM, mat);
ctVals = sampleVolumeTrilinearCountGuided(Img, vox);
blobVals = sampleVolumeTrilinearCountGuided(blob, vox);
ctScore = clipCountGuided((ctVals - p95) ./ max(p99_8 - p95, 1e-6), 0, 2.0);
blobScore = clipCountGuided(blobVals, 0, 2.0);
vals = 0.60 * ctScore + 0.40 * blobScore;
end

function [gridS, gridPts] = extendGridByEvidenceAutoCount(polylinePts, gridS, spacingMm, ...
    Img, blob, mat, p95, p99_8, ProjSurfRaw, supportDist, outsideDist, ...
    maxExtraTotal, baseEvThr, strongEvThr)
gridS = gridS(:);
gridPts = samplePolylineExtrapCountGuided(polylinePts, gridS);
for it = 1:maxExtraTotal
    leftS = gridS(1) - spacingMm;
    leftPt = samplePolylineExtrapCountGuided(polylinePts, leftS);
    leftEv = contactEvidenceValuesAutoCount(Img, blob, mat, p95, p99_8, leftPt);
    leftIn = LeG_intriangulation(ProjSurfRaw.vertices, ProjSurfRaw.faces, leftPt) == 1;
    leftDist = nearestSurfaceDistanceCountGuided(leftPt, ProjSurfRaw.vertices);
    leftOk = (leftEv >= baseEvThr) && (leftIn || leftDist <= supportDist || leftEv >= strongEvThr);
    leftBad = (~leftIn) && (leftDist > outsideDist) && (leftEv < strongEvThr);

    rightS = gridS(end) + spacingMm;
    rightPt = samplePolylineExtrapCountGuided(polylinePts, rightS);
    rightEv = contactEvidenceValuesAutoCount(Img, blob, mat, p95, p99_8, rightPt);
    rightIn = LeG_intriangulation(ProjSurfRaw.vertices, ProjSurfRaw.faces, rightPt) == 1;
    rightDist = nearestSurfaceDistanceCountGuided(rightPt, ProjSurfRaw.vertices);
    rightOk = (rightEv >= baseEvThr) && (rightIn || rightDist <= supportDist || rightEv >= strongEvThr);
    rightBad = (~rightIn) && (rightDist > outsideDist) && (rightEv < strongEvThr);

    optS = zeros(0, 1);
    optEv = zeros(0, 1);
    optPts = zeros(0, 3);
    optSide = zeros(0, 1);
    if leftOk && ~leftBad
        optS(end+1, 1) = leftS; %#ok<AGROW>
        optEv(end+1, 1) = leftEv; %#ok<AGROW>
        optPts(end+1, :) = leftPt; %#ok<AGROW>
        optSide(end+1, 1) = 1; %#ok<AGROW>
    end
    if rightOk && ~rightBad
        optS(end+1, 1) = rightS; %#ok<AGROW>
        optEv(end+1, 1) = rightEv; %#ok<AGROW>
        optPts(end+1, :) = rightPt; %#ok<AGROW>
        optSide(end+1, 1) = 2; %#ok<AGROW>
    end
    if isempty(optEv)
        break;
    end
    [~, bi] = max(optEv);
    if optSide(bi) == 1
        gridS = [optS(bi); gridS]; %#ok<AGROW>
        gridPts = [optPts(bi, :); gridPts]; %#ok<AGROW>
    else
        gridS = [gridS; optS(bi)]; %#ok<AGROW>
        gridPts = [gridPts; optPts(bi, :)]; %#ok<AGROW>
    end
end
end

function shafts = inflateSelectedShaftsAutoCount(shafts, Img, blob, mat, ProjSurfRaw, surfaceSupport)
if isempty(shafts)
    return;
end

v = Img(isfinite(Img));
pct = prctile(v, [95 99.8]);
p95 = pct(1);
p99_8 = pct(2);
[supportDist, outsideDist, surfaceWeight, baseEvThr, strongEvThr] = ...
    surfaceDistanceParamsAutoCount(surfaceSupport);
out = emptyAutoCountShaftStruct();

for i = 1:numel(shafts)
    shaft = shafts(i);
    if surfaceSupport >= 0.50
        cleanHighSupport = shaft.expectedContacts <= 12 && ...
            shaft.brainSupport >= 0.90 && shaft.outsideFrac <= 0.10;
        if ~cleanHighSupport
            out(end+1) = shaft; %#ok<AGROW>
            continue;
        end
    elseif surfaceSupport >= 0.42
        cleanMidSupport = shaft.expectedContacts <= 14 && ...
            shaft.brainSupport >= 0.80 && shaft.outsideFrac <= 0.18 && ...
            (shaft.supportScore >= 0.58 || shaft.evidenceScore >= 0.96);
        if ~cleanMidSupport
            out(end+1) = shaft; %#ok<AGROW>
            continue;
        end
    end

    maxExtraTotal = 0;
    if shaft.expectedContacts >= 5
        maxExtraTotal = 1;
        if shaft.evidenceScore >= 1.05
            maxExtraTotal = maxExtraTotal + 2;
        end
        if shaft.evidenceScore >= 1.25
            maxExtraTotal = maxExtraTotal + 1;
        end
        if shaft.supportScore >= 0.72
            maxExtraTotal = maxExtraTotal + 1;
        end
        if shaft.brainSupport >= 0.28
            maxExtraTotal = maxExtraTotal + 1;
        end
        if shaft.expectedContacts >= 8
            maxExtraTotal = maxExtraTotal + 1;
        end
        if surfaceSupport < 0.45
            maxExtraTotal = maxExtraTotal + 1;
        end
    end

    maxExtraCap = 8;
    if surfaceSupport >= 0.50
        maxExtraCap = 1;
    elseif surfaceSupport >= 0.42
        maxExtraCap = 2;
    elseif surfaceSupport >= 0.35
        maxExtraCap = 5;
    end
    maxExtraTotal = min(max(maxExtraTotal, 0), maxExtraCap);
    if maxExtraTotal == 0
        out(end+1) = shaft; %#ok<AGROW>
        continue;
    end

    gridS = (0:size(shaft.contactsMM, 1)-1)' * shaft.spacingMm;
    evRelax = 0;
    strongRelax = 0;
    if surfaceSupport < 0.42 && shaft.evidenceScore >= 1.15
        evRelax = evRelax + 0.06;
    end
    if surfaceSupport < 0.35 && shaft.supportScore >= 0.82
        evRelax = evRelax + 0.03;
    end
    if surfaceSupport < 0.35 && shaft.evidenceScore >= 1.30
        strongRelax = strongRelax + 0.06;
    end
    tunedBaseEv = max(0.50, baseEvThr - evRelax);
    tunedStrongEv = max(tunedBaseEv + 0.16, strongEvThr - strongRelax);
    [gridSNew, gridPtsNew] = extendGridByEvidenceAutoCount(shaft.contactsMM, gridS, ...
        shaft.spacingMm, Img, blob, mat, p95, p99_8, ProjSurfRaw, supportDist, ...
        outsideDist, maxExtraTotal, tunedBaseEv, tunedStrongEv);
    if size(gridPtsNew, 1) <= size(shaft.contactsMM, 1)
        out(end+1) = shaft; %#ok<AGROW>
        continue;
    end
    out(end+1) = rebuildShaftFromContactsAutoCount(shaft, gridPtsNew, Img, blob, ... %#ok<AGROW>
        mat, p95, p99_8, ProjSurfRaw, supportDist, outsideDist, surfaceWeight);
end

if isempty(out)
    shafts = out;
    return;
end
[~, ord] = sort([out.score], 'descend');
shafts = out(ord);
end

function shaft = rebuildShaftFromContactsAutoCount(shaftIn, contactsMM, Img, blob, mat, ...
    p95, p99_8, ProjSurfRaw, supportDist, outsideDist, surfaceWeight)
[~, ~, tSorted, orderedPts] = shaftAxisAutoCount(contactsMM);
orderedPts = deduplicateContactsAlongAxisAutoCount(orderedPts, tSorted, 1.5);
spacingMm = clipCountGuided(estimateSpacingCountGuided(orderedPts), 2.4, 5.0);
orderedPts = fillSmallContactGapsAutoCount(orderedPts, 1.6, 2.6);
evidenceVals = contactEvidenceValuesAutoCount(Img, blob, mat, p95, p99_8, orderedPts);
evidenceScore = mean(evidenceVals);
inside = LeG_intriangulation(ProjSurfRaw.vertices, ProjSurfRaw.faces, orderedPts) == 1;
dist = nearestSurfaceDistanceCountGuided(orderedPts, ProjSurfRaw.vertices);
brainSupport = mean(inside | (dist <= supportDist));
outsideFrac = mean(~inside & dist > outsideDist);
insideFrac = mean(inside);
regularity = spacingRegularityCountGuided(orderedPts, spacingMm);
addedContacts = max(0, size(orderedPts, 1) - size(shaftIn.contactsMM, 1));
supportScore = clipCountGuided(0.70 * shaftIn.supportScore + 0.30 * regularity, 0, 1);
score = shaftIn.score + 0.18 * addedContacts + ...
    0.30 * (evidenceScore - shaftIn.evidenceScore) + ...
    0.20 * (supportScore - shaftIn.supportScore) + ...
    surfaceWeight * (0.25 * (brainSupport - shaftIn.brainSupport) - ...
    0.35 * (outsideFrac - shaftIn.outsideFrac));

shaft = shaftIn;
shaft.expectedContacts = size(orderedPts, 1);
shaft.score = score;
shaft.supportScore = supportScore;
shaft.evidenceScore = evidenceScore;
shaft.spacingMm = spacingMm;
shaft.brainSupport = brainSupport;
shaft.insideFrac = insideFrac;
shaft.outsideFrac = outsideFrac;
shaft.contactsMM = orderedPts;
end

function shafts = rescueIntracranialFitsAutoCount(selected, fits, scoreThr, surfaceSupport)
if isempty(fits)
    shafts = selected;
    return;
end

shafts = selected;
rescue = emptyAutoCountShaftStruct();
targetScore = scoreThr - 0.70;
if surfaceSupport < 0.40
    targetScore = targetScore - 0.10;
end

for i = 1:numel(fits)
    shaft = fits(i);
    if shaft.expectedContacts < 4
        continue;
    end
    strongCore = shaft.insideFrac >= 0.22 || ...
        (shaft.brainSupport >= 0.46 && shaft.outsideFrac <= 0.26);
    if ~strongCore
        continue;
    end
    if shaft.outsideFrac > 0.42 && shaft.insideFrac < 0.18
        continue;
    end
    if shaft.supportScore < 0.24 || shaft.evidenceScore < 0.30
        continue;
    end
    if shaft.score < targetScore && shaft.evidenceScore < 0.75
        continue;
    end
    if isDuplicateShaftAutoCountGuided(shafts, shaft) || isDuplicateShaftAutoCountGuided(rescue, shaft)
        continue;
    end
    rescue(end+1) = shaft; %#ok<AGROW>
end

if ~isempty(rescue)
    rescueScore = [rescue.score]' + 1.00 * [rescue.insideFrac]' + ...
        0.45 * [rescue.brainSupport]' - 0.80 * [rescue.outsideFrac]';
    [~, ord] = sort(rescueScore, 'descend');
    rescue = rescue(ord);
    if isempty(shafts)
        shafts = rescue;
    else
        shafts = [shafts rescue]; %#ok<AGROW>
    end
end
end

function shafts = trimOutsideShaftEndsAutoCount(shafts, Img, blob, mat, ProjSurfRaw, surfaceSupport)
if isempty(shafts)
    return;
end

v = Img(isfinite(Img));
pct = prctile(v, [95 99.8]);
p95 = pct(1);
p99_8 = pct(2);
[supportDist, outsideDist, surfaceWeight, ~, strongEvThr] = ...
    surfaceDistanceParamsAutoCount(surfaceSupport);
trimDist = min(outsideDist, 5.5);
keepDist = min(outsideDist, trimDist + 1.5);
out = emptyAutoCountShaftStruct();
dropWholeOutside = 0;
trimmedAny = 0;

for i = 1:numel(shafts)
    shaft = shafts(i);
    pts = shaft.contactsMM;
    if size(pts, 1) < 3
        out(end+1) = shaft; %#ok<AGROW>
        continue;
    end

    evidenceVals = contactEvidenceValuesAutoCount(Img, blob, mat, p95, p99_8, pts);
    inside = LeG_intriangulation(ProjSurfRaw.vertices, ProjSurfRaw.faces, pts) == 1;
    dist = nearestSurfaceDistanceCountGuided(pts, ProjSurfRaw.vertices);
    insideCount = sum(inside);
    farOutside = ~inside & dist > (trimDist + 2.0);
    weakEv = evidenceVals < strongEvThr;
    if insideCount == 0 && size(pts, 1) >= 4 && mean(farOutside) >= 0.75 && mean(weakEv) >= 0.40
        dropWholeOutside = dropWholeOutside + 1;
        continue;
    end
    if insideCount <= 1 && size(pts, 1) >= 6 && mean(farOutside) >= 0.85 && mean(weakEv) >= 0.50
        dropWholeOutside = dropWholeOutside + 1;
        continue;
    end
    keep = inside | dist <= trimDist | ((dist <= keepDist) & (evidenceVals >= strongEvThr));
    keepIdx = find(keep);
    if isempty(keepIdx)
        continue;
    end

    lo = keepIdx(1);
    hi = keepIdx(end);
    if lo == 1 && hi == size(pts, 1)
        out(end+1) = shaft; %#ok<AGROW>
        continue;
    end

    leftTrim = lo - 1;
    rightTrim = size(pts, 1) - hi;
    if max(leftTrim, rightTrim) < 2 && shaft.brainSupport >= 0.80
        out(end+1) = shaft; %#ok<AGROW>
        continue;
    end

    if (hi - lo + 1) < 3
        continue;
    end

    shaftTrim = rebuildShaftFromContactsAutoCount(shaft, pts(lo:hi, :), Img, blob, ...
        mat, p95, p99_8, ProjSurfRaw, supportDist, outsideDist, surfaceWeight);
    if shaftTrim.insideFrac < 0.05 && size(shaftTrim.contactsMM, 1) >= 4
        dropWholeOutside = dropWholeOutside + 1;
        continue;
    end
    if shaftTrim.outsideFrac > 0.40 && shaftTrim.brainSupport < 0.30
        continue;
    end
    if lo > 1 || hi < size(pts, 1)
        trimmedAny = trimmedAny + 1;
    end
    out(end+1) = shaftTrim; %#ok<AGROW>
end

if isempty(out)
    shafts = out;
    return;
end
[~, ord] = sort([out.score], 'descend');
shafts = out(ord);
if dropWholeOutside > 0
    fprintf('[AutoElec] Outside reject: dropped %d shafts with no intracranial core\n', ...
        dropWholeOutside);
end
if trimmedAny > 0
    fprintf('[AutoElec] Outside end trim touched %d shafts\n', trimmedAny);
end
end

function [supportDist, outsideDist, surfaceWeight, baseEvThr, strongEvThr] = ...
    surfaceDistanceParamsAutoCount(surfaceSupport)
supportDist = 3.5;
outsideDist = 6.0;
surfaceWeight = 1.0;
baseEvThr = 0.64;
strongEvThr = 0.94;
if surfaceSupport < 0.25
    supportDist = 8.0;
    outsideDist = 11.0;
    surfaceWeight = 0.25;
    baseEvThr = 0.56;
    strongEvThr = 0.82;
elseif surfaceSupport < 0.35
    supportDist = 6.5;
    outsideDist = 9.5;
    surfaceWeight = 0.45;
    baseEvThr = 0.58;
    strongEvThr = 0.86;
elseif surfaceSupport < 0.45
    supportDist = 5.0;
    outsideDist = 8.0;
    surfaceWeight = 0.70;
    baseEvThr = 0.60;
    strongEvThr = 0.90;
end
end

function shafts = mergeCollinearShaftsAutoCount(shafts, strictMode)
if isempty(shafts)
    return;
end
changed = true;
while changed
    changed = false;
    for i = 1:numel(shafts)-1
        for j = i+1:numel(shafts)
            if ~shouldMergeShaftsAutoCount(shafts(i), shafts(j), strictMode)
                continue;
            end
            shafts(i) = combineShaftsAutoCount(shafts(i), shafts(j));
            shafts(j) = [];
            changed = true;
            break;
        end
        if changed
            break;
        end
    end
end
if ~isempty(shafts)
    [~, ord] = sort([shafts.score], 'descend');
    shafts = shafts(ord);
end
end

function shafts = dropRedundantFragmentsAutoCount(shafts)
if isempty(shafts)
    return;
end
[~, ord] = sortrows([-[shafts.expectedContacts]' -[shafts.score]']);
kept = emptyAutoCountShaftStruct();
for ii = 1:numel(ord)
    shaft = shafts(ord(ii));
    redundant = false;
    for j = 1:numel(kept)
        if isRedundantFragmentAutoCount(shaft, kept(j))
            redundant = true;
            break;
        end
    end
    if redundant
        continue;
    end
    kept(end+1) = shaft; %#ok<AGROW>
end
if ~isempty(kept)
    [~, ord] = sort([kept.score], 'descend');
    shafts = kept(ord);
else
    shafts = kept;
end
end

function shafts = consolidateShaftGroupsAutoCount(shafts, Img, blob, mat, ProjSurfRaw, surfaceSupport, strictMode)
if numel(shafts) < 2
    return;
end

nS = numel(shafts);
adj = false(nS, nS);
for i = 1:nS-1
    for j = i+1:nS
        if shouldGroupShaftsAutoCount(shafts(i), shafts(j), strictMode)
            adj(i, j) = true;
            adj(j, i) = true;
        end
    end
end

G = graph(adj);
comp = conncomp(G);
uComp = unique(comp(:))';
out = emptyAutoCountShaftStruct();

for k = 1:numel(uComp)
    members = find(comp == uComp(k));
    if numel(members) == 1
        out(end+1) = shafts(members); %#ok<AGROW>
        continue;
    end

    counts = [shafts(members).expectedContacts]';
    scores = [shafts(members).score]';
    [~, ord] = sortrows([-counts -scores]);
    orderedMembers = members(ord);

    if strictMode
        anchor = shafts(orderedMembers(1));
        for ii = 2:numel(orderedMembers)
            frag = shafts(orderedMembers(ii));
            if shouldAttachFragmentAutoCount(anchor, frag)
                anchor = combineShaftsAutoCount(anchor, frag);
            end
        end
        out(end+1) = anchor; %#ok<AGROW>
        continue;
    end

    fragmentedMidSupport = surfaceSupport >= 0.40 && surfaceSupport <= 0.52 && ...
        numel(shafts) >= 24 && max(counts) <= 10;
    if fragmentedMidSupport
        merged = shafts(orderedMembers(1));
        for ii = 2:numel(orderedMembers)
            merged = combineShaftsAutoCount(merged, shafts(orderedMembers(ii)));
        end
        out(end+1) = merged; %#ok<AGROW>
        continue;
    end

    pts = zeros(0, 3);
    seedIdx = zeros(1, 0);
    for ii = 1:numel(orderedMembers)
        pts = [pts; shafts(orderedMembers(ii)).contactsMM]; %#ok<AGROW>
        seedIdx = unique([seedIdx shafts(orderedMembers(ii)).seedIdx]); %#ok<AGROW>
    end
    chain = syntheticChainFromContactsAutoCount(pts, seedIdx, Img, blob, mat, ProjSurfRaw);
    if isempty(chain)
        merged = shafts(orderedMembers(1));
        for ii = 2:numel(orderedMembers)
            merged = combineShaftsAutoCount(merged, shafts(orderedMembers(ii)));
        end
        out(end+1) = merged; %#ok<AGROW>
        continue;
    end

    refit = fitContactsToChainAutoCount(chain, Img, blob, mat, ProjSurfRaw, surfaceSupport);
    if isempty(refit)
        merged = shafts(orderedMembers(1));
        for ii = 2:numel(orderedMembers)
            merged = combineShaftsAutoCount(merged, shafts(orderedMembers(ii)));
        end
        out(end+1) = merged; %#ok<AGROW>
    else
        out(end+1) = refit; %#ok<AGROW>
    end
end

if ~isempty(out)
    [~, ord] = sort([out.score], 'descend');
    shafts = out(ord);
else
    shafts = out;
end
end

function tf = shouldGroupShaftsAutoCount(a, b, strictMode)
[axisA, muA, ~, ptsA] = shaftAxisAutoCount(a.contactsMM);
[axisB, muB, ~, ptsB] = shaftAxisAutoCount(b.contactsMM);
angle = acosd(min(1, max(-1, abs(dot(axisA, axisB)))));
lineDist = shaftLineDistanceAutoCount(muA, axisA, muB, axisB);
if strictMode
    maxLineDist = 2.6;
    angleCap = 22.0;
    relaxedAngleCap = 26.0;
    relaxedLineDist = 1.8;
    maxGap = 10.0;
else
    maxLineDist = 4.5;
    angleCap = 28.0;
    relaxedAngleCap = 35.0;
    relaxedLineDist = 3.5;
    maxGap = 16.0;
end
if lineDist > maxLineDist
    tf = false;
    return;
end
if angle > angleCap
    tf = angle <= relaxedAngleCap && lineDist <= relaxedLineDist;
    if ~tf
        return;
    end
end

spacingClose = abs(a.spacingMm - b.spacingMm) <= 1.8 || ...
    min(a.expectedContacts, b.expectedContacts) <= 4;
if ~spacingClose
    tf = false;
    return;
end

gap = shaftIntervalGapAutoCount(ptsA, ptsB, axisA, axisB);
if gap > maxGap
    tf = false;
    return;
end

endsA = [ptsA(1, :); ptsA(end, :)];
endsB = [ptsB(1, :); ptsB(end, :)];
endDist = min(min(pdist2(endsA, endsB)));
smallPair = max(a.expectedContacts, b.expectedContacts) <= 6 && ...
    min(a.expectedContacts, b.expectedContacts) <= 4;
if smallPair && lineDist <= 1.5 && endDist <= 6.0 && gap <= 6.0
    tf = true;
    return;
end

spacingRef = clipCountGuided(median([a.spacingMm b.spacingMm]), 2.4, 5.0);
[meanNN, medNN] = shaftContactAlignmentAutoCount(a.contactsMM, b.contactsMM);
if gap <= 0.5 && min(a.expectedContacts, b.expectedContacts) >= 8
    if strictMode
        if medNN > max(2.2, 0.70 * spacingRef) || meanNN > max(4.0, 1.20 * spacingRef)
            tf = false;
            return;
        end
    else
        if medNN > max(5.5, 1.70 * spacingRef) || meanNN > max(10.0, 3.00 * spacingRef)
            tf = false;
            return;
        end
    end
end

if strictMode && min(a.expectedContacts, b.expectedContacts) >= 10
    tf = lineDist <= 1.8 && gap <= 6.0 && (endDist <= 6.0 || gap <= 4.5);
elseif strictMode
    tf = endDist <= 8.0 || gap <= 8.0;
else
    tf = endDist <= 14.0 || gap <= 14.0;
end
end

function tf = isRedundantFragmentAutoCount(fragment, other)
if fragment.expectedContacts > 5
    tf = false;
    return;
end
if other.expectedContacts < fragment.expectedContacts + 4
    tf = false;
    return;
end

[axisF, muF, ~, ptsF] = shaftAxisAutoCount(fragment.contactsMM);
[axisO, muO, ~, ptsO] = shaftAxisAutoCount(other.contactsMM);
angle = acosd(min(1, max(-1, abs(dot(axisF, axisO)))));
if angle > 14.0
    tf = false;
    return;
end

lineDist = shaftLineDistanceAutoCount(muF, axisF, muO, axisO);
if lineDist > 2.8
    tf = false;
    return;
end

gap = shaftIntervalGapAutoCount(ptsF, ptsO, axisF, axisO);
maxGap = max(7.5, 2.0 * max(fragment.spacingMm, other.spacingMm));
if gap > maxGap
    tf = false;
    return;
end

endF = [ptsF(1, :); ptsF(end, :)];
endO = [ptsO(1, :); ptsO(end, :)];
endDist = min(min(pdist2(endF, endO)));
if endDist > max(8.0, 2.4 * max(fragment.spacingMm, other.spacingMm)) && ...
        gap > max(5.0, 1.5 * max(fragment.spacingMm, other.spacingMm))
    tf = false;
    return;
end

D = pdist2(ptsF, ptsO);
meanMinContactDist = mean(min(D, [], 2));
tf = meanMinContactDist <= 3.2;
end

function chain = syntheticChainFromContactsAutoCount(pointsMM, seedIdx, Img, blob, mat, ProjSurfRaw)
chain = [];
if size(pointsMM, 1) < 3
    return;
end

v = Img(isfinite(Img));
pct = prctile(v, [95 99.8]);
p95 = pct(1);
p99_8 = pct(2);
[~, ~, tSorted, orderedPts] = shaftAxisAutoCount(pointsMM);
orderedPts = deduplicateContactsAlongAxisAutoCount(orderedPts, tSorted, 1.5);
if size(orderedPts, 1) < 3
    return;
end

spacing = estimateSpacingCountGuided(orderedPts);
regularity = spacingRegularityCountGuided(orderedPts, spacing);
span = polylineLengthCountGuided(orderedPts);
if regularity < 0.12 || span < 6.0
    return;
end

evidenceVals = contactEvidenceValuesAutoCount(Img, blob, mat, p95, p99_8, orderedPts);
inside = LeG_intriangulation(ProjSurfRaw.vertices, ProjSurfRaw.faces, orderedPts) == 1;
dist = nearestSurfaceDistanceCountGuided(orderedPts, ProjSurfRaw.vertices);
near = dist <= 8.0;
brainSupport = mean(inside | near);
insideFrac = mean(inside);
nearFrac = mean(near);
score = 1.4 * regularity + 0.75 * mean(evidenceVals) + ...
    0.20 * size(orderedPts, 1) + 0.05 * span + 0.35 * brainSupport;

chain.idx = seedIdx(:)';
chain.pointsMM = orderedPts;
chain.scores = evidenceVals(:);
chain.spacing = clipCountGuided(spacing, 2.4, 5.0);
chain.span = span;
chain.regularity = regularity;
chain.brainSupport = brainSupport;
chain.insideFrac = insideFrac;
chain.nearFrac = nearFrac;
chain.score = score;
end

function tf = shouldMergeShaftsAutoCount(a, b, strictMode)
[axisA, muA, ~, ptsA] = shaftAxisAutoCount(a.contactsMM);
[axisB, muB, ~, ptsB] = shaftAxisAutoCount(b.contactsMM);
angle = acosd(min(1, max(-1, abs(dot(axisA, axisB)))));
smallPair = max(a.expectedContacts, b.expectedContacts) <= 6 && ...
    min(a.expectedContacts, b.expectedContacts) <= 4;
if smallPair && angle <= 65.0
    lineDist = shaftLineDistanceAutoCount(muA, axisA, muB, axisB);
    if lineDist > 1.5
        tf = false;
        return;
    end
    gap = shaftIntervalGapAutoCount(ptsA, ptsB, axisA, axisB);
    if gap > 6.0
        tf = false;
        return;
    end
    endA = [ptsA(1, :); ptsA(end, :)];
    endB = [ptsB(1, :); ptsB(end, :)];
    endDist = min(min(pdist2(endA, endB)));
    tf = endDist <= 6.0;
    return;
end

if strictMode
    angleCap = 16.0;
    maxLineDist = 2.0;
    maxGap = 8.0;
else
    angleCap = 22.0;
    maxLineDist = 3.0;
    maxGap = 12.0;
end
if angle > angleCap
    tf = false;
    return;
end

lineDist = shaftLineDistanceAutoCount(muA, axisA, muB, axisB);
if lineDist > maxLineDist
    tf = false;
    return;
end

gap = shaftIntervalGapAutoCount(ptsA, ptsB, axisA, axisB);
if gap > maxGap
    tf = false;
    return;
end

endA = [ptsA(1, :); ptsA(end, :)];
endB = [ptsB(1, :); ptsB(end, :)];
endDist = min(min(pdist2(endA, endB)));
spacingRef = clipCountGuided(median([a.spacingMm b.spacingMm]), 2.4, 5.0);
[meanNN, medNN] = shaftContactAlignmentAutoCount(a.contactsMM, b.contactsMM);
overlapAxis = commonIntervalOverlapAutoCount(a.contactsMM, b.contactsMM, axisA, axisB);
if overlapAxis > 0.5
    if strictMode
        if medNN > max(2.2, 0.70 * spacingRef) || meanNN > max(4.0, 1.20 * spacingRef)
            tf = false;
            return;
        end
        tf = angle <= 8.0 && lineDist <= 1.2;
        return;
    end
    if medNN > max(5.5, 1.70 * spacingRef) || meanNN > max(10.0, 3.00 * spacingRef)
        tf = false;
        return;
    end
    tf = angle <= 16.0 && lineDist <= 2.2;
    return;
end

if strictMode && min(a.expectedContacts, b.expectedContacts) >= 10
    tf = lineDist <= 1.5 && gap <= 4.5 && endDist <= max(4.5, 1.35 * spacingRef);
elseif strictMode
    tf = endDist <= max(6.0, 1.7 * spacingRef);
else
    tf = endDist <= max(9.0, 2.1 * spacingRef);
end
end

function [meanNN, medNN] = shaftContactAlignmentAutoCount(pointsA, pointsB)
D = pdist2(pointsA, pointsB);
minA = min(D, [], 2);
minB = min(D, [], 1)';
combined = [minA; minB];
meanNN = mean(combined);
medNN = median(combined);
end

function overlap = commonIntervalOverlapAutoCount(pointsA, pointsB, axisA, axisB)
sgn = sign(dot(axisA, axisB));
if sgn == 0
    sgn = 1;
end
axisM = axisA + sgn * axisB;
axisM = axisM ./ max(norm(axisM), 1e-6);
tA = pointsA * axisM;
tB = pointsB * axisM;
overlap = min(max(tA), max(tB)) - max(min(tA), min(tB));
end

function tf = shouldAttachFragmentAutoCount(anchor, frag)
if frag.expectedContacts > max(7, anchor.expectedContacts - 4)
    tf = false;
    return;
end

[axisA, muA, ~, ~] = shaftAxisAutoCount(anchor.contactsMM);
[axisB, muB, ~, ~] = shaftAxisAutoCount(frag.contactsMM);
angle = acosd(min(1, max(-1, abs(dot(axisA, axisB)))));
if angle > 12.0
    tf = false;
    return;
end

lineDist = shaftLineDistanceAutoCount(muA, axisA, muB, axisB);
spacingRef = clipCountGuided(median([anchor.spacingMm frag.spacingMm]), 2.4, 5.0);
if lineDist > max(1.2, 0.35 * spacingRef)
    tf = false;
    return;
end

gap = shaftIntervalGapAutoCount(anchor.contactsMM, frag.contactsMM, axisA, axisB);
overlapAxis = commonIntervalOverlapAutoCount(anchor.contactsMM, frag.contactsMM, axisA, axisB);
[~, medNN] = shaftContactAlignmentAutoCount(anchor.contactsMM, frag.contactsMM);
endA = [anchor.contactsMM(1, :); anchor.contactsMM(end, :)];
endB = [frag.contactsMM(1, :); frag.contactsMM(end, :)];
endDist = min(min(pdist2(endA, endB)));

if overlapAxis > 0.5
    tf = false;
    return;
end
if gap > max(4.8, 1.55 * spacingRef)
    tf = false;
    return;
end
if endDist > max(5.5, 1.75 * spacingRef)
    tf = false;
    return;
end
if medNN > max(6.0, 1.9 * spacingRef)
    tf = false;
    return;
end
tf = true;
end

function shaft = combineShaftsAutoCount(a, b)
allPts = [a.contactsMM; b.contactsMM];
[~, ~, tSorted, orderedPts] = shaftAxisAutoCount(allPts);
mergeSpacing = clipCountGuided(median([a.spacingMm b.spacingMm]), 2.4, 5.0);
dedupSpacing = clipCountGuided(0.62 * mergeSpacing, 1.8, 2.4);
mergedPts = deduplicateContactsAlongAxisAutoCount(orderedPts, tSorted, dedupSpacing);
mergedPts = fillSmallContactGapsAutoCount(mergedPts, 1.8, 2.35);
spacingMm = estimateSpacingCountGuided(mergedPts);
score = max(a.score, b.score) + 0.18 * min(size(mergedPts, 1), 12) / 12.0;
supportScore = weightedAverageAutoCount([a.supportScore b.supportScore], ...
    [size(a.contactsMM,1) size(b.contactsMM,1)]);
evidenceScore = weightedAverageAutoCount([a.evidenceScore b.evidenceScore], ...
    [size(a.contactsMM,1) size(b.contactsMM,1)]);
brainSupport = weightedAverageAutoCount([a.brainSupport b.brainSupport], ...
    [size(a.contactsMM,1) size(b.contactsMM,1)]);
insideFrac = weightedAverageAutoCount([a.insideFrac b.insideFrac], ...
    [size(a.contactsMM,1) size(b.contactsMM,1)]);
outsideFrac = weightedAverageAutoCount([a.outsideFrac b.outsideFrac], ...
    [size(a.contactsMM,1) size(b.contactsMM,1)]);

shaft.expectedContacts = size(mergedPts, 1);
shaft.observedContacts = a.observedContacts + b.observedContacts;
shaft.score = score;
shaft.supportScore = supportScore;
shaft.evidenceScore = evidenceScore;
shaft.spacingMm = spacingMm;
shaft.brainSupport = brainSupport;
shaft.insideFrac = insideFrac;
shaft.outsideFrac = outsideFrac;
shaft.contactsMM = mergedPts;
shaft.seedIdx = unique([a.seedIdx b.seedIdx]);
end

function [axis1, mu, tSorted, orderedPts] = shaftAxisAutoCount(pointsMM)
mu = mean(pointsMM, 1);
if size(pointsMM, 1) < 2
    axis1 = [1; 0; 0];
else
    [~, ~, V] = svd(pointsMM - mu, 'econ');
    axis1 = V(:, 1);
end
t = (pointsMM - mu) * axis1;
[tSorted, ord] = sort(t);
orderedPts = pointsMM(ord, :);
end

function d = shaftLineDistanceAutoCount(muA, axisA, muB, axisB)
cr = cross(axisA(:), axisB(:));
crNorm = norm(cr);
if crNorm < 1e-6
    d = norm(cross((muB - muA)', axisA(:)));
else
    d = abs(dot((muB - muA), cr')) / crNorm;
end
end

function gap = shaftIntervalGapAutoCount(pointsA, pointsB, axisA, axisB)
sgn = sign(dot(axisA(:), axisB(:)));
if sgn == 0
    sgn = 1;
end
axisM = axisA(:) + sgn * axisB(:);
axisM = axisM ./ max(norm(axisM), 1e-6);
tA = pointsA * axisM;
tB = pointsB * axisM;
gap = max(0, max(min(tA) - max(tB), min(tB) - max(tA)));
end

function pointsMM = deduplicateContactsAlongAxisAutoCount(pointsMM, tSorted, minSpacingMm)
if isempty(pointsMM)
    return;
end
keepPts = pointsMM(1, :);
keepT = tSorted(1);
for i = 2:size(pointsMM, 1)
    if tSorted(i) - keepT(end) < minSpacingMm
        keepPts(end, :) = 0.5 * (keepPts(end, :) + pointsMM(i, :));
        keepT(end, 1) = 0.5 * (keepT(end) + tSorted(i));
    else
        keepPts(end+1, :) = pointsMM(i, :); %#ok<AGROW>
        keepT(end+1, 1) = tSorted(i); %#ok<AGROW>
    end
end
pointsMM = keepPts;
end

function pointsMM = fillSmallContactGapsAutoCount(pointsMM, gapLoFactor, gapHiFactor)
if size(pointsMM, 1) < 3
    return;
end
[~, ~, tSorted, orderedPts] = shaftAxisAutoCount(pointsMM);
gaps = diff(tSorted);
valid = gaps(gaps > 2.0 & gaps < 6.0);
if isempty(valid)
    pointsMM = orderedPts;
    return;
end
spacingMm = median(valid);
out = orderedPts(1, :);
for i = 1:numel(gaps)
    gap = gaps(i);
    if gap > gapLoFactor * spacingMm && gap < gapHiFactor * spacingMm
        nMissing = max(0, round(gap / spacingMm) - 1);
        for m = 1:nMissing
            alpha = m / (nMissing + 1);
            out(end+1, :) = (1 - alpha) * orderedPts(i, :) + alpha * orderedPts(i+1, :); %#ok<AGROW>
        end
    end
    out(end+1, :) = orderedPts(i+1, :); %#ok<AGROW>
end
pointsMM = out;
end

function val = weightedAverageAutoCount(values, weights)
weights = double(weights(:));
values = double(values(:));
if sum(weights) <= 0
    val = mean(values);
else
    val = sum(values .* weights) ./ sum(weights);
end
end

function shafts = emptyAutoCountShaftStruct()
shafts = struct('expectedContacts', {}, 'observedContacts', {}, ...
    'score', {}, 'supportScore', {}, 'spacingMm', {}, ...
    'contactsMM', {}, 'seedIdx', {}, 'brainSupport', {}, ...
    'insideFrac', {}, 'evidenceScore', {}, 'outsideFrac', {});
end

function shafts = selectAndFitShaftsCountGuided(chains, expectedCounts, Img, blob, mat)
shafts = struct('expectedContacts', {}, 'observedContacts', {}, ...
    'score', {}, 'supportScore', {}, 'spacingMm', {}, ...
    'contactsMM', {}, 'seedIdx', {}, 'brainSupport', {}, ...
    'insideFrac', {}, 'evidenceScore', {}, 'outsideFrac', {});
if isempty(chains)
    return;
end

nExp = numel(expectedCounts);
nChains = min(numel(chains), max(nExp * 5, 40));
chains = chains(1:nChains);
cost = buildAssignmentCostCountGuided(chains, expectedCounts);

usedRows = false(nChains, 1);
usedCols = false(nExp, 1);
usedMask = false(1, max([chains.idx]));
selectedCounts = zeros(0, 1);

[~, ord] = sort(cost(:), 'ascend');
for ii = 1:numel(ord)
    [r, c] = ind2sub(size(cost), ord(ii));
    if r > nChains || c > nExp || usedRows(r) || usedCols(c)
        continue;
    end
    overlap = mean(usedMask(chains(r).idx));
    if overlap > 0.72
        continue;
    end
    shaft = fitContactsToChainCountGuided(chains(r), expectedCounts(c), Img, blob, mat);
    if isempty(shaft)
        continue;
    end
    shafts(end+1) = shaft; %#ok<AGROW>
    usedRows(r) = true;
    usedCols(c) = true;
    usedMask(chains(r).idx) = true;
    selectedCounts(end+1,1) = expectedCounts(c); %#ok<AGROW>
    if numel(shafts) == nExp
        break;
    end
end

remaining = expectedCounts(:);
for i = 1:numel(selectedCounts)
    k = find(remaining == selectedCounts(i), 1, 'first');
    if ~isempty(k)
        remaining(k) = [];
    end
end

for i = 1:numel(remaining)
    expCount = remaining(i);
    bestRow = 0;
    bestCost = Inf;
    for r = 1:nChains
        if usedRows(r)
            continue;
        end
        overlap = mean(usedMask(chains(r).idx));
        if overlap > 0.82
            continue;
        end
        cst = singleAssignmentCostCountGuided(chains(r), expCount);
        if cst < bestCost
            bestCost = cst;
            bestRow = r;
        end
    end
    if bestRow == 0
        continue;
    end
    shaft = fitContactsToChainCountGuided(chains(bestRow), expCount, Img, blob, mat);
    if isempty(shaft)
        continue;
    end
    shafts(end+1) = shaft; %#ok<AGROW>
    usedRows(bestRow) = true;
    usedMask(chains(bestRow).idx) = true;
end

if ~isempty(shafts)
    [~, ord] = sort([shafts.score], 'descend');
    shafts = shafts(ord);
    shafts = shafts(1:min(nExp, numel(shafts)));
end
end

function cost = buildAssignmentCostCountGuided(chains, expectedCounts)
n = max(numel(chains), numel(expectedCounts));
cost = 40 * ones(n, n);
for i = 1:numel(chains)
    for j = 1:numel(expectedCounts)
        cost(i, j) = singleAssignmentCostCountGuided(chains(i), expectedCounts(j));
    end
end
end

function cst = singleAssignmentCostCountGuided(chain, expCount)
baseCount = max(2, round(chain.span ./ max(chain.spacing, 1e-3)) + 1);
obsCount = numel(chain.idx);
mismatch = abs(baseCount - expCount);
obsGap = max(0, expCount - obsCount);
cst = 1.4 * mismatch + 0.55 * obsGap ...
    - 1.15 * chain.regularity - 0.08 * chain.span ...
    - 0.14 * obsCount - 0.12 * chain.score;
end

function shaft = fitContactsToChainCountGuided(chain, expectedCount, Img, blob, mat)
shaft = [];
pts = chain.pointsMM;
obsScores = chain.scores;
obsS = arcLengthsCountGuided(pts);
spacingGuess = min(max(chain.spacing, 2.6), 4.6);
totalLen = max(obsS(end), spacingGuess);

v = Img(isfinite(Img));
pct = prctile(v, [95 99.8]);
p95 = pct(1);
p99_8 = pct(2);

bestScore = -Inf;
bestSupport = 0;
bestSpacing = spacingGuess;
bestPts = zeros(0, 3);

spacingGrid = max(2.5, spacingGuess - 0.7):0.15:min(4.8, spacingGuess + 0.7);
for spacing = spacingGrid
    nominalLen = spacing * max(expectedCount - 1, 1);
    startMin = min(obsS(1), totalLen - nominalLen) - 4.0;
    startMax = max(obsS(end) - nominalLen, 0.0) + 4.0;
    for s0 = startMin:0.5:startMax
        gridS = s0 + spacing * (0:expectedCount-1)';
        gridPts = samplePolylineExtrapCountGuided(pts, gridS);
        vox = mmToVoxCountGuided(gridPts, mat);
        ctVals = sampleVolumeTrilinearCountGuided(Img, vox);
        blobVals = sampleVolumeTrilinearCountGuided(blob, vox);
        ctScore = clipCountGuided((ctVals - p95) ./ max(p99_8 - p95, 1e-6), 0, 2.0);
        blobScore = clipCountGuided(blobVals, 0, 2.0);

        dObs = abs(bsxfun(@minus, gridS(:), obsS(:)'));
        support = mean(exp(-((min(dObs, [], 1) ./ 1.1).^2)));
        evidence = mean(0.65 * ctScore + 0.35 * blobScore);
        score = 2.0 * support + 1.4 * evidence + 1.2 * chain.regularity ...
            + 0.15 * mean(obsScores) - 0.08 * abs(spacing - spacingGuess);
        if score > bestScore
            bestScore = score;
            bestSupport = support;
            bestSpacing = spacing;
            bestPts = gridPts;
        end
    end
end

if isempty(bestPts)
    return;
end

shaft.expectedContacts = expectedCount;
shaft.observedContacts = numel(chain.idx);
shaft.score = bestScore;
shaft.supportScore = bestSupport;
shaft.spacingMm = bestSpacing;
shaft.contactsMM = bestPts;
shaft.seedIdx = chain.idx;
shaft.brainSupport = chain.brainSupport;
shaft.insideFrac = chain.insideFrac;
shaft.evidenceScore = mean(contactEvidenceValuesAutoCount(Img, blob, mat, p95, p99_8, bestPts));
shaft.outsideFrac = max(0, 1 - chain.brainSupport);
end

function s = arcLengthsCountGuided(pointsMM)
s = zeros(size(pointsMM, 1), 1);
if size(pointsMM, 1) > 1
    s(2:end) = cumsum(sqrt(sum(diff(pointsMM, 1, 1).^2, 2)));
end
end

function pts = samplePolylineExtrapCountGuided(pointsMM, sQuery)
obsS = arcLengthsCountGuided(pointsMM);
pts = zeros(numel(sQuery), 3);
if size(pointsMM, 1) == 1
    pts = repmat(pointsMM, numel(sQuery), 1);
    return;
end

startDir = pointsMM(2, :) - pointsMM(1, :);
startDir = startDir ./ max(norm(startDir), 1e-6);
endDir = pointsMM(end, :) - pointsMM(end-1, :);
endDir = endDir ./ max(norm(endDir), 1e-6);

for i = 1:numel(sQuery)
    sq = sQuery(i);
    if sq <= 0
        pts(i, :) = pointsMM(1, :) + startDir .* sq;
        continue;
    end
    if sq >= obsS(end)
        pts(i, :) = pointsMM(end, :) + endDir .* (sq - obsS(end));
        continue;
    end
    idx = find(obsS <= sq, 1, 'last');
    idx = min(max(idx, 1), size(pointsMM, 1) - 1);
    segLen = max(obsS(idx+1) - obsS(idx), 1e-6);
    alpha = (sq - obsS(idx)) ./ segLen;
    pts(i, :) = (1 - alpha) .* pointsMM(idx, :) + alpha .* pointsMM(idx+1, :);
end
end

function mm = voxToMMCountGuided(vox, mat)
mmH = [vox, ones(size(vox,1),1)] * mat';
mm = mmH(:, 1:3);
end

function vox = mmToVoxCountGuided(mm, mat)
voxH = mat \ [mm, ones(size(mm,1),1)]';
vox = voxH(1:3, :)';
end

function vals = sampleVolumeTrilinearCountGuided(V, vox)
if isempty(vox)
    vals = zeros(0, 1);
    return;
end
vox = double(vox);
vox(:, 1) = min(max(vox(:, 1), 1), size(V, 1));
vox(:, 2) = min(max(vox(:, 2), 1), size(V, 2));
vox(:, 3) = min(max(vox(:, 3), 1), size(V, 3));
vals = interpn(V, vox(:,1), vox(:,2), vox(:,3), 'linear', 0);
vals = vals(:);
end

function x = clipCountGuided(x, lo, hi)
x(x < lo) = lo;
x(x > hi) = hi;
end

function diagPlotCountGuided(app, ThrHU, NumObj, idx, elapsedTime)
try
    fH = figure('Position', [50 50 400 400], 'Name', app.PatientIDStr);
    aH = axes('Parent', fH);
    plot(aH, ThrHU, NumObj, '-');
    hold(aH, 'on');
    if idx >= 1 && idx <= numel(ThrHU)
        plot(aH, ThrHU(idx), NumObj(idx), 'or');
    end
    xlabel(aH, 'Threshold (HU)');
    ylabel(aH, '# Detections');
    title(aH, sprintf('%0.1f s chain-fit', elapsedTime));
catch
end
end
