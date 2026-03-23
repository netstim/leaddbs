function [WC, T] = LeG_autoElecs(app, varargin)
StartTime = tic;

fprintf('\n[AutoElec] ========== Starting electrode detection ==========\n');

% === Unpack inputs ========================================================
Img         = single(app.CTImg);
MaxElecs    = getDetectOpt(app, 'maxelecs', 250, [10, 5000]);
MRInfo      = app.MRInfo;
ProjSurfRaw = app.ProjSurfRaw;
imSz        = size(Img);

% Voxel dimensions in mm  [dx  dy  dz]
voxDim = sqrt(sum(app.CTInfo.mat(1:3,1:3).^2));

fprintf('[AutoElec] Image size: %d x %d x %d, voxel dims: %.2f x %.2f x %.2f mm\n', ...
    imSz(1), imSz(2), imSz(3), voxDim(1), voxDim(2), voxDim(3));
fprintf('[AutoElec] MaxElecs = %d\n', MaxElecs);
det = getDetectStruct(app);

% Phase flags: which pipeline stages to run (default all on for backward compat)
if isfield(det, 'usePhase2')
    usePhase2 = logical(det.usePhase2);
else
    usePhase2 = true;
end
if isfield(det, 'usePhase3')
    usePhase3 = logical(det.usePhase3);
else
    usePhase3 = true;
end
if isfield(det, 'usePhase4')
    usePhase4 = logical(det.usePhase4);
else
    usePhase4 = false;
end

% Raw threshold sweep kept for the backward-compatible diagnostic plot
TMax  = (app.CTRng(4)-app.CTRng(1)) / (app.CTRng(2)-app.CTRng(1));
TMin  = 1;
ThrR  = linspace(TMax, TMin, 21);  ThrR(1) = [];
ThrHU = ThrR*(app.CTRng(2)-app.CTRng(1)) + app.CTRng(1);
ThrHU = (ThrHU - app.CTInfo.pinfo(2)) ./ app.CTInfo.pinfo(1);

% Optional detector mode:
%   'chainfit' (default): deterministic count-free peak/chain detector.
%       Bright compact CT peaks seed shaft hypotheses; the surface mesh is
%       used only for candidate gating and shaft/contact scoring. This is
%       the primary maintained path.
%   'chainoptimization': experimental global shaft-selection variant.
%       Uses the same candidate/fit stages as chainfit, but selects shafts
%       via a beam-search conflict optimizer before refinement.
%   'countguided' (explicit only): same fitter with expected per-electrode
%       contact counts as priors, when such priors are intentionally supplied.
%   'chainhybrid': legacy rescue fusion kept only as an explicit fallback.
%   'hybrid': fuse shaft-first + legacy blob detections.
%   'shaftfirst': tube/shaft-first pipeline + contact fitting.
%   'legacy': existing blob+RANSAC pipeline below.
detectorMode = 'chainoptimization';
if isfield(det, 'method') && ~isempty(det.method)
    if isstring(det.method) || ischar(det.method)
        detectorMode = lower(char(det.method));
    end
elseif isfield(det, 'mode') && ~isempty(det.mode)
    if isstring(det.mode) || ischar(det.mode)
        detectorMode = lower(char(det.mode));
    end
end
if nargin >= 2 && ~isempty(varargin{1})
    fm = varargin{1};
    if ischar(fm) || isstring(fm)
        detectorMode = lower(char(fm));
    end
end

if strcmp(detectorMode, 'chainhybrid')
    [WC, T] = runChainHybridPipeline(app, Img, MRInfo, ProjSurfRaw, ...
                                     voxDim, ThrR, ThrHU, TMax, TMin, StartTime);
    return;
elseif strcmp(detectorMode, 'countguided') || strcmp(detectorMode, 'chainfit')
    [WC, T] = LeG_autoElecs_countguided(app, Img, MRInfo, ProjSurfRaw, ...
        voxDim, ThrR, ThrHU, TMax, TMin, StartTime, detectorMode);
    return;
elseif strcmp(detectorMode, 'chainoptimization')
    [WC, T] = LeG_autoElecs_chainoptimization(app, Img, MRInfo, ProjSurfRaw, ...
        voxDim, ThrR, ThrHU, TMax, TMin, StartTime);
    return;
end

if strcmp(detectorMode, 'hybrid')
    fprintf('[AutoElec] Detector mode: hybrid\n');
    [WC, T] = runHybridPipeline(app, Img, MRInfo, ProjSurfRaw, ...
                                voxDim, ThrR, ThrHU, TMax, TMin, StartTime);
    return;
elseif strcmp(detectorMode, 'shaftfirst')
    fprintf('[AutoElec] Detector mode: shaftfirst\n');
    [WC, T] = runShaftFirstPipeline(app, Img, MRInfo, ProjSurfRaw, ...
                                    voxDim, ThrR, ThrHU, TMax, TMin, StartTime);
    return;
end
fprintf('[AutoElec] Detector mode: legacy\n');

% =========================================================================
%  PHASE 1 — MULTI-SCALE BLOB DETECTION
% =========================================================================

% 1a  Robust intensity percentiles (ignore background <= 0) ---------------
vValid = Img(Img > 0);
if isempty(vValid)
    WC = []; T = 0;
    diagPlot(app, ThrHU, zeros(size(ThrHU)), 1, toc(StartTime));
    return;
end
pct  = prctile(vValid(:), [50 75 90 95 99 99.9]);
pMed = pct(1);  p75 = pct(2);  p90 = pct(3);  p95 = pct(4);  p99 = pct(5);  p99_9 = pct(6);

fprintf('[AutoElec] Intensity percentiles — p50=%.1f  p75=%.1f  p90=%.1f  p95=%.1f  p99=%.1f  p99.9=%.1f\n', ...
    pct(1), pct(2), pct(3), pct(4), pct(5), pct(6));

% 1b  Laplacian-of-Gaussian via Difference-of-Gaussians -------------------
%     DoG ≈ sigma^2 * LoG.  Anisotropic smoothing (sigma per axis) so that
%     thick-slice / oblique electrodes still get a strong blob response.
%     Sigma in mm, then converted to voxel units per axis.
nScales = 4;
sigmas_mm = linspace(0.4, 1.2, nScales);

blob = zeros(imSz, 'single');
for s = 1:nScales
    sigma_mm = sigmas_mm(s);
    sigma_vox = sigma_mm ./ voxDim(:)';           % [sx sy sz] in voxel units per axis
    g1 = imgaussfilt3(Img, sigma_vox);
    g2 = imgaussfilt3(Img, sigma_vox * 1.6);
    blob = max(blob, -(g2 - g1) * sigma_mm^2);
end

bClip = prctile(blob(blob > 0), 99.9);
if bClip > 0
    blob = blob / bClip;
end

fprintf('[AutoElec] LoG blob detection: %d scales, anisotropic (sigma_mm = %.2f – %.2f)\n', ...
    nScales, sigmas_mm(1), sigmas_mm(end));
fprintf('[AutoElec]   Blob response range: [%.3f, %.3f], clip = %.3f\n', ...
    min(blob(:)), max(blob(:)), bClip);

% 1c  Combined score: normalised intensity x blob response -----------------
%     Scale by p99.9 so the very brightest voxels can exceed 1 (no cap).
%     That avoids missing the brightest contacts when LoG is weak (metal bloom).
%     The floor of 0.1 on blob still helps weak-blob bright spots.
intNorm   = max(0, (Img - pMed) / max(p99_9 - pMed, 1e-6));
combScore = intNorm .* max(blob, 0.1);

% 1d  Volume bounds for a single contact -----------------------------------
%     sEEG contact geometry: cylinder ~0.8 mm diameter, ~1.5 mm long
%     Physical volume ≈ pi * 0.4^2 * 1.5 ≈ 0.75 mm^3
%     Generous bounds accommodate partial-volume effects and metal bloom.
voxVol = prod(voxDim);
minVox = max(2, round(0.75 * 0.15 / voxVol));
maxVox = max(60, round(0.75 * 45 / voxVol));

fprintf('[AutoElec] Volume bounds: [%d, %d] voxels (voxVol = %.3f mm^3)\n', ...
    minVox, maxVox, voxVol);

% 1e  Threshold sweep on the combined score --------------------------------
%     By sweeping on the combined intensity-blob score (rather than raw
%     intensity alone), detection is largely independent of absolute HU
%     calibration.  Candidates from ALL thresholds are collected and merged.
nThr    = 18;
thrLevs = linspace(0.05, 0.9, nThr);

% Accumulator: each row = [row  col  slice  meanIntensity  blobVal  thrIdx]
allCand = zeros(0, 6);

fprintf('[AutoElec] Sweeping %d thresholds on combined score [%.2f – %.2f]...\n', ...
    nThr, thrLevs(1), thrLevs(end));

for k = 1:nThr
    CC = bwconncomp(combScore > thrLevs(k), 26);
    if CC.NumObjects == 0, continue; end

    nRawCC = CC.NumObjects;
    ccSz = cellfun(@numel, CC.PixelIdxList);
    ccSz = ccSz(:);   % ensure column for element-wise logic
    % Get MeanIntensity for all CCs so we can keep very bright small blobs
    propsAll = regionprops3(CC, Img, ...
        'WeightedCentroid', 'MeanIntensity', 'PrincipalAxisLength');
    miAll = propsAll.MeanIntensity(:);   % ensure column, same length as ccSz
    keep = (ccSz >= minVox & ccSz <= maxVox) | ...
           (ccSz >= 1 & ccSz <= 3 & miAll >= p99_9);
    CC.PixelIdxList = CC.PixelIdxList(keep);
    CC.NumObjects   = sum(keep);
    nVolOK = CC.NumObjects;
    if CC.NumObjects == 0, continue; end

    % Shape filter — reject very elongated blobs (streak artefacts)
    pal = propsAll.PrincipalAxisLength(keep, :);     % [N x 3] descending
    if size(pal, 2) == 3
        elongation = pal(:,1) ./ max(pal(:,3), 0.1);
        goodShape  = elongation < 8;
    else
        goodShape  = true(CC.NumObjects, 1);
    end
    nShapeOK = sum(goodShape);

    wc = propsAll.WeightedCentroid(keep, :);
    wc = wc(goodShape, :);                           % [col  row  slice]
    mi = propsAll.MeanIntensity(keep);
    mi = mi(goodShape);
    if isempty(wc), continue; end

    wc(:, [1 2]) = wc(:, [2 1]);                     % → [row  col  slice]

    % Keep only detections inside the brain / skull envelope, or very bright
    % (likely electrode at boundary — don't drop brightest contacts)
    wcMM = [wc, ones(size(wc,1),1)] * MRInfo.mat';
    wcMM(:,4) = [];
    inBrain = LeG_intriangulation(ProjSurfRaw.vertices, ...
                                  ProjSurfRaw.faces, wcMM);
    keepBrain = inBrain | (mi >= p99_9);
    wc = wc(keepBrain, :);
    mi = mi(keepBrain);
    nBrainOK = size(wc, 1);
    if isempty(wc), continue; end

    % Look up the blob value at each centroid
    rc = max(min(round(wc), imSz), 1);
    linIdx = sub2ind(imSz, rc(:,1), rc(:,2), rc(:,3));
    bv = blob(linIdx);

    allCand = [allCand; ...
        wc, double(mi), double(bv), k*ones(size(wc,1),1)]; %#ok<AGROW>

    fprintf('[AutoElec]   Thr %2d (%.2f): %d raw CC -> %d vol-OK -> %d shape-OK -> %d in-brain\n', ...
        k, thrLevs(k), nRawCC, nVolOK, nShapeOK, nBrainOK);
end

fprintf('[AutoElec] Total raw candidates across all thresholds: %d\n', size(allCand, 1));

% --- Raw-threshold detection count curve (for the diagnostic plot) --------
NumObj = zeros(numel(ThrR), 1);
for kt = 1:numel(ThrR)
    cc0  = bwconncomp(Img > ThrR(kt), 26);
    cs0  = cellfun(@numel, cc0.PixelIdxList);
    NumObj(kt) = sum(cs0 >= minVox & cs0 <= maxVox);
end

if isempty(allCand)
    WC = []; T = 0;
    diagPlot(app, ThrHU, NumObj, round(numel(ThrR)/2), toc(StartTime));
    return;
end

% 1f  Merge duplicate detections across thresholds -------------------------
%     The same physical contact is typically detected at many thresholds.
%     Agglomerative clustering with a 1.5 mm cut-off merges these, and the
%     final centroid is a confidence-weighted average of the group.
candVox  = allCand(:, 1:3);
candMM   = [candVox, ones(size(candVox,1),1)] * MRInfo.mat';
candMM(:,4) = [];
candInt  = allCand(:, 4);
candBlob = allCand(:, 5);

% Pre-filter: keep top 3000 candidates (by intensity*blob), but never drop
% the very brightest (intensity >= p99.9) so we don't miss bright electrodes.
maxNCand = 3000;
nBeforeCap = size(candMM, 1);
if size(candMM, 1) > maxNCand
    protected = (candInt >= p99_9);
    nProt = sum(protected);
    rawScore = candInt .* max(candBlob, 0.01);
    if nProt >= maxNCand
        selProt = find(protected);
        [~, sIdx] = sort(rawScore(protected), 'descend');
        sel = selProt(sIdx(1:maxNCand));
    else
        selProt = find(protected);
        rest = find(~protected);
        [~, sIdxRest] = sort(rawScore(rest), 'descend');
        nTake = min(maxNCand - nProt, numel(rest));
        sel = [selProt; rest(sIdxRest(1:nTake))];
    end
    candMM  = candMM(sel, :);
    candVox = candVox(sel, :);
    candInt = candInt(sel);
    candBlob = candBlob(sel);
    allCand  = allCand(sel, :);
    fprintf('[AutoElec] Pre-filter: capped %d candidates to %d (top by score; %d brightest protected)\n', ...
        nBeforeCap, size(candMM, 1), nProt);
end

if size(candMM, 1) > 1
    Z   = linkage(candMM, 'average');
    cID = cluster(Z, 'cutoff', 1.5, 'criterion', 'distance');
else
    cID = 1;
end

nGroups = max(cID);
mrgMM   = zeros(nGroups, 3);
mrgVox  = zeros(nGroups, 3);
mrgInt  = zeros(nGroups, 1);
mrgBlob = zeros(nGroups, 1);
mrgCnt  = zeros(nGroups, 1);

for g = 1:nGroups
    mask = cID == g;
    w = candInt(mask) .* max(candBlob(mask), 0.01);
    w = w / sum(w);
    mrgMM(g, :)  = w' * candMM(mask, :);
    mrgVox(g, :) = w' * candVox(mask, :);
    mrgInt(g)    = max(candInt(mask));
    mrgBlob(g)   = max(candBlob(mask));
    mrgCnt(g)    = sum(mask);
end

% Confidence score: persistence across thresholds + intensity + blob shape
detectFrac = mrgCnt / nThr;
intScore   = (mrgInt - min(mrgInt)) / max(max(mrgInt) - min(mrgInt), 1e-6);
confidence = 0.4 * detectFrac + 0.3 * intScore + 0.3 * mrgBlob;

fprintf('[AutoElec] Cross-threshold merge: %d candidates -> %d unique contacts\n', ...
    nBeforeCap, nGroups);
fprintf('[AutoElec] Confidence distribution: min=%.2f  median=%.2f  max=%.2f  (mean=%.2f)\n', ...
    min(confidence), median(confidence), max(confidence), mean(confidence));
fprintf('[AutoElec]   Contacts with conf > 0.5: %d,  conf > 0.3: %d,  conf <= 0.3: %d\n', ...
    sum(confidence > 0.5), sum(confidence > 0.3 & confidence <= 0.5), sum(confidence <= 0.3));
fprintf('[AutoElec] Phase 1 complete (%.1f s)\n', toc(StartTime));

% ---------- Optional early exit: Phase 1 only (merged blob contacts) ------
% Optionally keep only contacts with confidence > phase1ConfMin (e.g. 0.5) at this stage.
if ~usePhase2
    fprintf('\n[AutoElec] Phase 2 disabled — returning merged contacts from Phase 1 only\n');
    phase1ConfMin = 0.5;   % no filter by default; set app.detect.phase1ConfMin = 0.5 to keep only conf > 0.5
    if isfield(app.detect, 'phase1ConfMin') && ~isempty(app.detect.phase1ConfMin)
        phase1ConfMin = double(app.detect.phase1ConfMin(1));
    end
    if phase1ConfMin > 0
        keepConf = (confidence > phase1ConfMin);
        nBefore = numel(confidence);
        mrgMM = mrgMM(keepConf, :);
        confidence = confidence(keepConf);
        fprintf('[AutoElec] Phase 1 confidence filter (conf > %.2f): %d -> %d contacts\n', ...
            phase1ConfMin, nBefore, numel(confidence));
    end
    if isempty(mrgMM)
        WC = []; T = (TMax + TMin) / 2;
        diagPlot(app, ThrHU, NumObj, round(numel(ThrR)/2), toc(StartTime));
        fprintf('\n[AutoElec] No contacts left after confidence filter. Exiting.\n\n');
        return;
    end
    [~, sortIdx] = sort(confidence, 'descend');
    finalMM_early = mrgMM(sortIdx, :);
    if size(finalMM_early, 1) > MaxElecs
        finalMM_early = finalMM_early(1:MaxElecs, :);
        fprintf('[AutoElec] Capped at MaxElecs: %d contacts\n', MaxElecs);
    end
    voxH = MRInfo.mat \ [finalMM_early, ones(size(finalMM_early, 1), 1)]';
    WC   = round(voxH(1:3, :)');
    T    = (TMax + TMin) / 2;
    EndTime = toc(StartTime);
    [~, idxD] = min(abs(NumObj - size(WC, 1)));
    diagPlot(app, ThrHU, NumObj, idxD, EndTime);
    if phase1ConfMin > 0
        fprintf('\n[AutoElec] ========== FINAL (Phase 1 only, conf > %.2f): %d contacts in %.1f s ==========\n\n', ...
            phase1ConfMin, size(WC, 1), EndTime);
    else
        fprintf('\n[AutoElec] ========== FINAL (Phase 1 only): %d contacts in %.1f s ==========\n\n', ...
            size(WC, 1), EndTime);
    end
    return;
end

% =========================================================================
%  PHASE 2 — RANSAC ELECTRODE-SHAFT DISCOVERY
% =========================================================================
%  sEEG electrodes are approximately straight: contacts are collinear.
%  Sequential RANSAC discovers these linear structures.  Each pair of
%  candidate points proposes a line; inliers within a 2 mm tube vote.
%  After two passes of PCA refinement the final inlier set is locked in
%  and removed before finding the next shaft.

nPts = size(mrgMM, 1);
if nPts < 3
    WC = round(mrgVox);
    T  = (TMax + TMin) / 2;
    diagPlot(app, ThrHU, NumObj, round(numel(ThrR)/2), toc(StartTime));
    return;
end

inlierRadius = 2.5;            % mm – tube radius for line inlier test (slight loosen for localization jitter)
minShaftPts  = 2;               % allow 2-point shafts so underdetected electrodes are kept
minSpanMm    = 6;               % require 2-point shafts to span at least this (avoid junk shafts)
maxShafts    = ceil(MaxElecs/3);
assigned     = false(nPts, 1);
shafts       = struct('idx',{}, 'dir',{}, 'mu',{});

rngOld = rng;
rng(42, 'twister');             % reproducible RANSAC sampling

fprintf('\n[AutoElec] === PHASE 2: RANSAC shaft discovery ===\n');
fprintf('[AutoElec] %d candidates to cluster, inlier radius = %.1f mm, min %d pts/shaft (2-pt span >= %.0f mm)\n', ...
    nPts, inlierRadius, minShaftPts, minSpanMm);

for sIter = 1:maxShafts
    uIdx = find(~assigned);
    nU   = numel(uIdx);
    if nU < minShaftPts, break; end

    ptsU  = mrgMM(uIdx, :);
    confU = confidence(uIdx);

    bestInliers = [];
    bestScore   = 0;
    nIter = min(5000, nU*(nU-1)/2);

    for it = 1:nIter
        % Bias sampling toward higher-confidence points so weaker shafts get found
        w = confU(:) + 1e-3;
        w = w / sum(w);
        i1 = randsample(nU, 1, true, w);
        i2 = randsample(nU, 1, true, w);
        if i1 == i2, continue; end
        pair = [i1, i2];
        lineDir = ptsU(pair(2),:) - ptsU(pair(1),:);
        lineLen = norm(lineDir);
        if lineLen < 2, continue; end              % too close — skip
        lineDir = lineDir / lineLen;

        vecs = ptsU - ptsU(pair(1), :);
        proj = vecs * lineDir';
        perpDist = sqrt(max(sum(vecs.^2, 2) - proj.^2, 0));

        inMask = perpDist < inlierRadius;
        if sum(inMask) < minShaftPts, continue; end
        % 2-point shafts: require minimum span so we don't create nonsense shafts
        if sum(inMask) == 2
            twoPts = ptsU(inMask, :);
            span = norm(twoPts(2,:) - twoPts(1,:));
            if span < minSpanMm, continue; end
        end

        score = sum(confU(inMask));
        if score > bestScore
            bestScore   = score;
            bestInliers = uIdx(inMask);
        end
    end

    if isempty(bestInliers), break; end

    % Two-pass PCA refinement: re-fit line and re-evaluate inliers
    for pass_ = 1:2
        sPts = mrgMM(bestInliers, :);
        mu   = mean(sPts, 1);
        [~, ~, V] = svd(sPts - mu, 'econ');
        lineDir = V(:,1)';

        vecsAll  = mrgMM(uIdx, :) - mu;
        projAll  = vecsAll * lineDir';
        perpAll  = sqrt(max(sum(vecsAll.^2, 2) - projAll.^2, 0));
        refined  = uIdx(perpAll < inlierRadius);
        if numel(refined) >= minShaftPts
            bestInliers = refined;
        end
    end

    % Final PCA on the locked-in inliers
    sPts = mrgMM(bestInliers, :);
    mu   = mean(sPts, 1);
    [~, ~, V] = svd(sPts - mu, 'econ');

    shafts(end+1) = struct('idx', bestInliers, ...
                           'dir', V(:,1)', ...
                           'mu',  mu);              %#ok<AGROW>
    assigned(bestInliers) = true;

    % Shaft diagnostics
    t_diag = (sPts - mu) * V(:,1);
    shaftLen = max(t_diag) - min(t_diag);
    meanConf = mean(confidence(bestInliers));
    fprintf('[AutoElec]   Shaft %d: %d contacts, span = %.1f mm, mean conf = %.2f, dir = [%.2f, %.2f, %.2f]\n', ...
        numel(shafts), numel(bestInliers), shaftLen, meanConf, V(1,1), V(2,1), V(3,1));
end

% Extend shafts: assign unassigned points that lie on an existing shaft line
% (within inlierRadius and near the shaft span); keep margin tight to avoid sparse false contacts
extendRadius = inlierRadius;
extendMargin = 5;   % mm beyond shaft endpoint to still count as "on shaft"
uList = find(~assigned);
for i = 1:numel(uList)
    u = uList(i);
    pt = mrgMM(u, :);
    bestS = [];
    bestPerp = inf;
    for s = 1:numel(shafts)
        sPts = mrgMM(shafts(s).idx, :);
        mu = mean(sPts, 1);
        dir = shafts(s).dir;
        tShaft = (sPts - mu) * dir';
        tMin = min(tShaft); tMax = max(tShaft);
        vec = pt - mu;
        proj = vec * dir';
        perpDist = sqrt(max(sum(vec.^2) - proj^2, 0));
        if perpDist < extendRadius && proj >= tMin - extendMargin && proj <= tMax + extendMargin
            if perpDist < bestPerp
                bestPerp = perpDist;
                bestS = s;
            end
        end
    end
    if ~isempty(bestS)
        shafts(bestS).idx = [shafts(bestS).idx; u];
        assigned(u) = true;
    end
end

rng(rngOld);                    % restore RNG state

nUnassigned = sum(~assigned);
fprintf('[AutoElec] RANSAC complete: %d shafts found, %d contacts assigned, %d unassigned (shaft extension applied)\n', ...
    numel(shafts), sum(assigned), nUnassigned);
fprintf('[AutoElec] Phase 2 complete (%.1f s)\n', toc(StartTime));

% =========================================================================
%  PHASE 3 — REGULAR-SPACING MODEL FIT & GAP FILLING (or shaft-only output)
% =========================================================================

finalMM     = zeros(0, 3);
finalConf   = zeros(0, 1);
finalShaftID = zeros(0, 1);   % 1..numel(shafts) per shaft, 0 = singleton

if ~usePhase3
    % Phase 3 disabled: use only RANSAC-assigned contacts per shaft (no grid/gap fill)
    fprintf('\n[AutoElec] Phase 3 disabled — using RANSAC shaft contacts only (no gap filling)\n');
    for s = 1:numel(shafts)
        sIdx = shafts(s).idx;
        nC = numel(sIdx);
        finalMM      = [finalMM;      mrgMM(sIdx, :)];  %#ok<AGROW>
        finalConf    = [finalConf;    confidence(sIdx)]; %#ok<AGROW>
        finalShaftID = [finalShaftID; s*ones(nC, 1)];  %#ok<AGROW>
    end
    unassigned = find(~assigned);
    nSingletons = 0;
    for u = 1:numel(unassigned)
        idxU = unassigned(u);
        if confidence(idxU) > 0.6 || mrgInt(idxU) >= p99_9
            finalMM      = [finalMM;      mrgMM(idxU, :)];       %#ok<AGROW>
            finalConf    = [finalConf;    confidence(idxU)];      %#ok<AGROW>
            finalShaftID = [finalShaftID; 0];                     %#ok<AGROW>
            nSingletons = nSingletons + 1;
        end
    end
    fprintf('[AutoElec] Shafts-only: %d contacts from shafts + %d singletons = %d total\n', ...
        size(finalMM, 1) - nSingletons, nSingletons, size(finalMM, 1));
    fprintf('[AutoElec] Phase 3 skipped (%.1f s)\n', toc(StartTime));
else
%  Phase 3 enabled: GAP-BASED filling only.
%  Purpose: We have N detected contacts. Real electrodes have regular
%  spacing; sometimes one contact is missed (weak blob). So we should ONLY
%  add interpolated contacts inside gaps that are clearly 2x or 3x the
%  spacing (one or two missed contacts). We must NOT lay a full grid over
%  the whole span — that wrongly assumes a contact every baseSp mm and
%  creates hundreds of false positives.
%  Steps:
%    3a. Estimate inter-contact spacing from observed gaps.
%    3b. For each consecutive pair of detected contacts, if gap ≈ k*baseSp
%        with k>=2, add (k-1) interpolated positions in that gap, each
%        validated by local intensity. Cap fill per gap and per shaft.

fprintf('\n[AutoElec] === PHASE 3: Gap-based fill only (no full-span grid) ===\n');

for s = 1:numel(shafts)
    sIdx = shafts(s).idx;
    sDir = shafts(s).dir;
    sMu  = shafts(s).mu;

    sPts  = mrgMM(sIdx, :);
    sConf = confidence(sIdx);

    % Project onto shaft axis
    t = (sPts - sMu) * sDir';
    [t, order] = sort(t);
    sPts  = sPts(order, :);
    sConf = sConf(order);
    nS    = numel(t);

    if nS < 2
        finalMM      = [finalMM;      sPts];   %#ok<AGROW>
        finalConf    = [finalConf;    sConf];   %#ok<AGROW>
        finalShaftID = [finalShaftID; s*ones(nS, 1)]; %#ok<AGROW>
        continue;
    end

    % --- 3a  Robust spacing estimation ------------------------------------
    gaps    = diff(t);
    posGaps = gaps(gaps > 0.5);          % ignore tiny residual gaps
    if isempty(posGaps)
        finalMM      = [finalMM;      sPts];   %#ok<AGROW>
        finalConf    = [finalConf;    sConf];   %#ok<AGROW>
        finalShaftID = [finalShaftID; s*ones(nS, 1)]; %#ok<AGROW>
        continue;
    end

    % Candidate spacings: observed gaps and their integer sub-divisions
    candSp = posGaps(:);
    for div = 2:4
        candSp = [candSp; posGaps(:)/div]; %#ok<AGROW>
    end
    candSp = candSp(candSp >= 1.5 & candSp <= 12);  % sEEG range (mm)

    if isempty(candSp)
        baseSp = median(posGaps);
    else
        scores = zeros(size(candSp));
        for cs = 1:numel(candSp)
            ratios   = posGaps / candSp(cs);
            rounded  = max(round(ratios), 1);
            residual = abs(ratios - rounded);
            scores(cs) = sum(residual < 0.25) - 0.5 * sum(residual);
        end
        [~, bestIdx] = max(scores);
        baseSp = candSp(bestIdx);
    end
    if baseSp < 1, baseSp = 3.5; end    % safe fallback

    % --- 3b  Gap-based fill only ------------------------------------------
    % Output list: start with first detected contact, then for each
    % consecutive pair (t_i, t_{i+1}), if gap suggests missed contacts
    % (gap in [1.5*baseSp, 2.5*baseSp] -> 1 missed, [2.5*baseSp, 3.5*baseSp] -> 2, etc.),
    % insert that many interpolated positions. Cap at maxFillPerGap
    % and ensure total contacts per shaft <= nS + reasonable.
    maxFillPerGap = 3;                  % allow up to 3 filled contacts per gap (dimmer mid-lead)
    maxTotalShaft  = nS + 8;           % total contacts on shaft at most nS+8

    shaftMM   = sPts(1, :);             % first contact
    shaftConf = sConf(1);

    nInterp   = 0;
    nInterpFail = 0;

    for i = 1:(nS - 1)
        gap = t(i+1) - t(i);
        % How many contacts are missing between t(i) and t(i+1)?
        % gap ≈ (nMissing+1)*baseSp  =>  nMissing = round(gap/baseSp) - 1
        nMissing = round(gap / baseSp) - 1;
        % Only fill if gap is clearly larger than one spacing (residual small)
        ratio = gap / baseSp;
        if nMissing < 1 || ratio < 1.35
            % No fill: gap is ~1 spacing or less
            shaftMM   = [shaftMM;   sPts(i+1, :)]; %#ok<AGROW>
            shaftConf = [shaftConf; sConf(i+1)];   %#ok<AGROW>
            continue;
        end
        nMissing = min(nMissing, maxFillPerGap);
        % Cap so shaft total stays <= maxTotalShaft (already have shaftMM; will add nMissing + rest of sPts)
        remainingDetected = nS - i;  % sPts(i+1)..sPts(nS) still to add
        if size(shaftMM, 1) + nMissing + remainingDetected > maxTotalShaft
            nMissing = max(0, maxTotalShaft - size(shaftMM, 1) - remainingDetected);
        end
        if nMissing < 1
            shaftMM   = [shaftMM;   sPts(i+1, :)]; %#ok<AGROW>
            shaftConf = [shaftConf; sConf(i+1)];   %#ok<AGROW>
            continue;
        end
        % Insert nMissing positions evenly between t(i) and t(i+1)
        for k = 1:nMissing
            tau = t(i) + (gap * k / (nMissing + 1));
            interpPt  = sMu + tau * sDir;
            interpVox = MRInfo.mat \ [interpPt'; 1];
            interpVox = round(interpVox(1:3))';

            if all(interpVox >= 1) && ...
               interpVox(1) <= imSz(1) && ...
               interpVox(2) <= imSz(2) && ...
               interpVox(3) <= imSz(3)

                [di, dj, dk] = ndgrid(-1:1, -1:1, -1:1);
                nbr = interpVox + [di(:), dj(:), dk(:)];
                nbr = max(nbr, 1);
                nbr(:,1) = min(nbr(:,1), imSz(1));
                nbr(:,2) = min(nbr(:,2), imSz(2));
                nbr(:,3) = min(nbr(:,3), imSz(3));
                nbrIdx    = sub2ind(imSz, nbr(:,1), nbr(:,2), nbr(:,3));
                localMean = mean(Img(nbrIdx));
                localBlob = mean(blob(nbrIdx));

                % Require clearly elevated intensity and blob; no dimmer tier to avoid spurious contacts
                if localMean > p90 && localBlob > 0.03
                    shaftMM   = [shaftMM;   interpPt]; %#ok<AGROW>
                    shaftConf = [shaftConf; 0.25];     %#ok<AGROW>
                    nInterp = nInterp + 1;
                else
                    nInterpFail = nInterpFail + 1;
                end
            else
                nInterpFail = nInterpFail + 1;
            end
        end

        shaftMM   = [shaftMM;   sPts(i+1, :)]; %#ok<AGROW>
        shaftConf = [shaftConf; sConf(i+1)];   %#ok<AGROW>
    end

    nDetected = nS;
    nTotal = size(shaftMM, 1);
    fprintf('[AutoElec]   Shaft %d: %d detected, %d gap-filled -> %d total (%d interp rejected)\n', ...
        s, nDetected, nInterp, nTotal, nInterpFail);

    finalMM      = [finalMM;      shaftMM];   %#ok<AGROW>
    finalConf    = [finalConf;    shaftConf];  %#ok<AGROW>
    finalShaftID = [finalShaftID; s*ones(nTotal, 1)]; %#ok<AGROW>
end

fprintf('[AutoElec] All shafts processed: %d contacts from shafts\n', size(finalMM, 1));

% Add singletons: only confident or very bright unassigned (avoid sparse false positives)
unassigned = find(~assigned);
nSingletons = 0;
for u = 1:numel(unassigned)
    idxU = unassigned(u);
    if confidence(idxU) > 0.6 || mrgInt(idxU) >= p99_9
        finalMM      = [finalMM;      mrgMM(idxU, :)];       %#ok<AGROW>
        finalConf    = [finalConf;    confidence(idxU)];      %#ok<AGROW>
        finalShaftID = [finalShaftID; 0];                     %#ok<AGROW>
        nSingletons = nSingletons + 1;
    end
end

fprintf('[AutoElec] Singletons: %d unassigned, %d added (conf > 0.6 or intensity >= p99.9)\n', ...
    numel(unassigned), nSingletons);
fprintf('[AutoElec] Entering Phase 4 with %d contacts total\n', size(finalMM, 1));
fprintf('[AutoElec] Phase 3 complete (%.1f s)\n', toc(StartTime));

end  % usePhase3

% Snapshot for diagnostics: contacts after Phase 3 (before bone/linearity filter)
finalMM_afterPhase3   = finalMM;
finalConf_afterPhase3 = finalConf;
finalShaftID_afterPhase3 = finalShaftID;

% =========================================================================
%  BONE / SURFACE FILTER — remove contacts outside brain and curved "shafts"
% =========================================================================
%  (1) Contact-level: load brain surface; drop every contact that lies
%      OUTSIDE the brain (bone, skull). Then drop shafts with < 2 contacts.
%  (2) Linearity: drop shafts whose contacts form a curve/arc (e.g. bone
%      following skull) rather than a straight line (real sEEG).
if size(finalMM, 1) > 0
    % --- (1) Contact-level brain filter ---
    surPath = app.SurfacesFile;
    loaded = load(surPath);
    if isfield(loaded, 'BrainSurfRaw')
        BrainSurfRaw = loaded.BrainSurfRaw;
    elseif isfield(loaded, 'k') && isfield(loaded.k, 'BrainSurfRaw')
        BrainSurfRaw = loaded.k.BrainSurfRaw;
    else
        fns = fieldnames(loaded);
        if ~isempty(fns) && isfield(loaded.(fns{1}), 'BrainSurfRaw')
            BrainSurfRaw = loaded.(fns{1}).BrainSurfRaw;
        else
            BrainSurfRaw = [];
        end
    end
    if ~isempty(BrainSurfRaw) && isfield(BrainSurfRaw, 'vertices') && isfield(BrainSurfRaw, 'faces')
        BrainV = BrainSurfRaw.vertices;
        BrainF = BrainSurfRaw.faces;
        % Per-contact: keep only contacts INSIDE the brain (drop bone)
        inBrain = LeG_intriangulation(BrainV, BrainF, finalMM);
        keep = inBrain;
        % Drop shafts that have fewer than 2 contacts left after removing bone
        shaftIDs = unique(finalShaftID(finalShaftID > 0));
        for s = shaftIDs(:)'
            if sum(keep & (finalShaftID == s)) < 2
                keep(finalShaftID == s) = false;
            end
        end
        % Singletons outside brain are also dropped (likely bone)
        nRemoved = sum(~keep);
        if nRemoved > 0
            finalMM      = finalMM(keep, :);
            finalConf    = finalConf(keep);
            finalShaftID = finalShaftID(keep);
            fprintf('[AutoElec] Brain filter: removed %d contacts outside brain (bone/surface)\n', nRemoved);
        end
    end
elseif isfield(app, 'SurfacesFile') && ~isempty(app.SurfacesFile)
    fprintf('[AutoElec] SurfacesFile not found or invalid, skipping brain filter\n');

    % --- (2) Linearity filter: drop shafts that form a curve/arc (bone) ---
    % Real sEEG contacts lie on a straight line; bone often follows skull contour.
    shaftIDs = unique(finalShaftID(finalShaftID > 0));
    if ~isempty(shaftIDs)
        maxCurvature = 0.25;   % max allowed (max perpendicular dist / span); above = curved → remove
        removeShaftCurve = false(max(shaftIDs), 1);
        for s = shaftIDs(:)'
            idx = (finalShaftID == s);
            pts = finalMM(idx, :);
            if size(pts, 1) < 3
                continue;
            end
            mu = mean(pts, 1);
            [~, ~, V] = svd(pts - mu, 'econ');
            proj = (pts - mu) * V(:, 1);
            span = max(proj) - min(proj);
            if span < 8
                continue;   % very short shaft, skip linearity check
            end
            perpDist = sqrt(max(sum((pts - mu).^2, 2) - proj.^2, 0));
            curvature = max(perpDist) / span;
            if curvature > maxCurvature
                removeShaftCurve(s) = true;
            end
        end
        keep = true(size(finalMM, 1), 1);
        for s = shaftIDs(:)'
            if removeShaftCurve(s)
                keep(finalShaftID == s) = false;
            end
        end
        nCurve = sum(~keep);
        if nCurve > 0
            finalMM      = finalMM(keep, :);
            finalConf    = finalConf(keep);
            finalShaftID = finalShaftID(keep);
            fprintf('[AutoElec] Linearity filter: removed %d contacts (%d curved shaft(s), likely bone)\n', ...
                nCurve, sum(removeShaftCurve));
        end
    end
end

% =========================================================================
%  PHASE 4 — DEDUPLICATION, RANKING, OUTPUT
% =========================================================================

fprintf('\n[AutoElec] === PHASE 4: Deduplication & output ===\n');

if isempty(finalMM)
    WC = []; T = 0;
    diagPlot(app, ThrHU, NumObj, round(numel(ThrR)/2), toc(StartTime));
    return;
end

nAfterBoneFilter = size(finalMM, 1);   % for diagnostic funnel (after bone/linearity filter)

% Sort by confidence (highest first) — always applied
[~, sortIdx] = sort(finalConf, 'descend');
finalMM      = finalMM(sortIdx, :);
finalConf    = finalConf(sortIdx);
finalShaftID = finalShaftID(sortIdx);

if usePhase4
    % Enforce minimum inter-contact distance (sEEG spacing is typically 3–5 mm;
    % avoid unrealistically close contacts by removing the lower-confidence of any pair too close)
    minInterContactMm = 2.0;
    nBeforeDedup = size(finalMM, 1);
    if size(finalMM, 1) > 1
        D = pdist2(finalMM, finalMM);
        D(logical(eye(size(D)))) = Inf;
        toRemove = false(size(finalMM, 1), 1);

        for i = 1:size(D, 1)
            if toRemove(i), continue; end
            tooClose = find(D(i, :) < minInterContactMm & ~toRemove');
            for j = tooClose
                if j <= i, continue; end
                if finalConf(j) <= finalConf(i)
                    toRemove(j) = true;
                else
                    toRemove(i) = true;
                    break;
                end
            end
        end

        nDedupRemoved = sum(toRemove);
        finalMM      = finalMM(~toRemove, :);
        finalConf    = finalConf(~toRemove);
        finalShaftID = finalShaftID(~toRemove);
        fprintf('[AutoElec] Min inter-contact %.1f mm: removed %d too-close contacts (%d -> %d)\n', ...
            minInterContactMm, nDedupRemoved, nBeforeDedup, size(finalMM, 1));
    end

    % Cap at MaxElecs
    nBeforeCap2 = size(finalMM, 1);
    if size(finalMM, 1) > MaxElecs
        finalMM      = finalMM(1:MaxElecs, :);
        finalConf    = finalConf(1:MaxElecs);
        finalShaftID = finalShaftID(1:MaxElecs);
        fprintf('[AutoElec] MaxElecs cap: trimmed %d -> %d contacts\n', ...
            nBeforeCap2, MaxElecs);
    end
else
    fprintf('[AutoElec] Phase 4 disabled — no dedup, no MaxElecs cap\n');
end

% Convert mm back to 1-based voxel coordinates
voxH = MRInfo.mat \ [finalMM, ones(size(finalMM, 1), 1)]';
WC   = round(voxH(1:3, :)');

% Nominal threshold for backward compatibility
T = (TMax + TMin) / 2;

% === Diagnostic plots ======================================================
EndTime = toc(StartTime);
[~, idxD] = min(abs(NumObj - size(WC, 1)));
diagPlot(app, ThrHU, NumObj, idxD, EndTime);
diagPlotsAutoElecs(app, ThrHU, NumObj, mrgMM, confidence, assigned, ...
    finalMM_afterPhase3, finalConf_afterPhase3, finalShaftID_afterPhase3, ...
    nAfterBoneFilter, finalMM, finalConf, finalShaftID, EndTime);

fprintf('\n[AutoElec] ========== FINAL: %d contacts detected in %.1f s ==========\n\n', ...
    size(WC, 1), EndTime);

end


% =========================================================================
%  Chain-hybrid pipeline (count-free chainfit + legacy blob rescue)
% =========================================================================
function [WC, T] = runChainHybridPipeline(app, Img, MRInfo, ProjSurfRaw, ...
                                          voxDim, ThrR, ThrHU, TMax, TMin, StartTime)
fprintf('\n[AutoElec-ChainHybrid] === Chain-fit + legacy fusion ===\n');

% Pass A: count-free chain-fit.
try
    [WCc, Tc] = LeG_autoElecs_countguided(app, Img, MRInfo, ProjSurfRaw, ...
        voxDim, ThrR, ThrHU, TMax, TMin, StartTime, 'chainfit');
    fprintf('[AutoElec-ChainHybrid] Chain-fit produced %d contacts\n', size(WCc,1));
catch
    warning('[AutoElec-ChainHybrid] Chain-fit pass failed, falling back to legacy+hybrid logic: %s');
    [WC, T] = runHybridPipeline(app, Img, MRInfo, ProjSurfRaw, ...
                                voxDim, ThrR, ThrHU, TMax, TMin, StartTime);
    return;
end

% Pass B: legacy blob+RANSAC.
try
    [WCl, Tl] = LeG_autoElecs(app, 'legacy');
    fprintf('[AutoElec-ChainHybrid] Legacy produced %d contacts\n', size(WCl,1));
catch
    warning('[AutoElec-ChainHybrid] Legacy pass failed: %s');
    WCl = zeros(0,3);
    Tl = Tc;
end

[WC, conf] = fuseHybridDetections(app, Img, MRInfo, ProjSurfRaw, voxDim, ...
                                  WCc, WCl, 'chainfit', 'legacy');
T = mean([Tc, Tl]);

NumObj = zeros(numel(ThrR), 1);
for kt = 1:numel(ThrR)
    cc0 = bwconncomp(Img > ThrR(kt), 26);
    if cc0.NumObjects == 0, continue; end
    cs0 = cellfun(@numel, cc0.PixelIdxList);
    NumObj(kt) = sum(cs0 >= 2);
end
[~, idxD] = min(abs(NumObj - size(WC,1)));
diagPlot(app, ThrHU, NumObj, idxD, toc(StartTime));

if isempty(conf), confMean = 0; else, confMean = mean(conf); end
fprintf('[AutoElec-ChainHybrid] FINAL fused contacts: %d (mean conf %.2f) in %.1f s\n', ...
    size(WC,1), confMean, toc(StartTime));
end


% =========================================================================
%  Hybrid pipeline (fuse shaft-first and legacy blob detections)
% =========================================================================
function [WC, T] = runHybridPipeline(app, Img, MRInfo, ProjSurfRaw, ...
                                     voxDim, ThrR, ThrHU, TMax, TMin, StartTime)
fprintf('\n[AutoElec-Hybrid] === Hybrid fusion (shaftfirst + legacy) ===\n');

% Pass A: shaft-first
[WCs, Ts] = runShaftFirstPipeline(app, Img, MRInfo, ProjSurfRaw, ...
                                  voxDim, ThrR, ThrHU, TMax, TMin, StartTime);
fprintf('[AutoElec-Hybrid] Shaft-first produced %d contacts\n', size(WCs,1));

% Pass B: legacy blob+RANSAC
try
    [WCl, Tl] = LeG_autoElecs(app, 'legacy');
    fprintf('[AutoElec-Hybrid] Legacy produced %d contacts\n', size(WCl,1));
catch ME
%     warning('[AutoElec-Hybrid] Legacy pass failed: %s', ME.message);
    WCl = zeros(0,3);
    Tl = Ts;
end

[WC, conf] = fuseHybridDetections(app, Img, MRInfo, ProjSurfRaw, voxDim, ...
                                  WCs, WCl, 'shaft', 'legacy');
T = mean([Ts, Tl]);

% Keep backward-compatible diagnostic curve display.
NumObj = zeros(numel(ThrR), 1);
for kt = 1:numel(ThrR)
    cc0 = bwconncomp(Img > ThrR(kt), 26);
    if cc0.NumObjects == 0, continue; end
    cs0 = cellfun(@numel, cc0.PixelIdxList);
    NumObj(kt) = sum(cs0 >= 2);
end
[~, idxD] = min(abs(NumObj - size(WC,1)));
diagPlot(app, ThrHU, NumObj, idxD, toc(StartTime));

if isempty(conf), confMean = 0; else, confMean = mean(conf); end
fprintf('[AutoElec-Hybrid] FINAL fused contacts: %d (mean conf %.2f) in %.1f s\n', ...
    size(WC,1), confMean, toc(StartTime));
end


% =========================================================================
%  Hybrid fusion helpers
% =========================================================================
function [WC, outConf] = fuseHybridDetections(app, Img, MRInfo, ProjSurfRaw, voxDim, WCs, WCl, primaryLabel, secondaryLabel)
imSz = size(Img);
MaxElecs = getDetectOpt(app, 'maxelecs', 250, [10, 5000]);
if nargin < 8 || isempty(primaryLabel), primaryLabel = 'primary'; end
if nargin < 9 || isempty(secondaryLabel), secondaryLabel = 'legacy'; end

if isempty(WCs) && isempty(WCl)
    WC = zeros(0,3);
    outConf = zeros(0,1);
    return;
end

candVox = [WCs; WCl];
src = [ones(size(WCs,1),1); 2*ones(size(WCl,1),1)]; % 1=shaft,2=legacy
inVol = candVox(:,1)>=1 & candVox(:,1)<=imSz(1) & ...
        candVox(:,2)>=1 & candVox(:,2)<=imSz(2) & ...
        candVox(:,3)>=1 & candVox(:,3)<=imSz(3);
candVox = candVox(inVol, :);
src = src(inVol);
if isempty(candVox)
    WC = zeros(0,3);
    outConf = zeros(0,1);
    return;
end

blob = computeBlobScoreMap(Img, voxDim);
[locI, locB] = sampleLocalAtVox(candVox, Img, blob);
vValid = Img(Img > 0);
pct = prctile(vValid(:), [90 99 99.8]);
p90 = pct(1); p99 = pct(2); p99_8 = pct(3);

iScore = min(max((locI - p90) / max(p99 - p90, 1e-6), 0), 1.5);
bPos = blob(blob > 0);
if isempty(bPos)
    bScale = 1;
else
    bScale = max(prctile(bPos, 95), 1e-6);
end
bScore = min(max(locB / bScale, 0), 1.5);
srcScore = (src == 1) * 0.64 + (src == 2) * 0.56;

candMM = voxToMM(candVox, MRInfo.mat);
[BrainV, BrainF] = loadBrainSurfaceForAutoElecs(app, ProjSurfRaw);
hasBrainSurf = ~isempty(BrainV) && ~isempty(BrainF);
nShaftIn = sum(src == 1);
nLegIn = sum(src == 2);
fprintf('[AutoElec-Hybrid] Fusion input: %s=%d %s=%d\n', primaryLabel, nShaftIn, secondaryLabel, nLegIn);

nearShaft = false(size(src));
if any(src == 1)
    dS = pdist2(candMM, candMM(src == 1, :));
    nearShaft = min(dS, [], 2) <= 4.2;
end

legacyRescue = false(size(src));
legacyChainScore = zeros(size(src));
if any(src == 2)
    legIdx = find(src == 2);
    [keepLeg, scoreLeg] = legacyChainSupport(candMM(legIdx, :), BrainV, BrainF, hasBrainSurf);
    legacyRescue(legIdx) = keepLeg;
    legacyChainScore(legIdx) = scoreLeg;
end

support = ones(size(src));
isLegacy = (src == 2);
support(isLegacy) = 0.44 + 0.34 * double(nearShaft(isLegacy)) + 0.30 * legacyChainScore(isLegacy);
weakLegacy = isLegacy & ~nearShaft & ~legacyRescue;
support(weakLegacy) = 0.22 + 0.22 * legacyChainScore(weakLegacy);
support = min(max(support, 0.15), 1.25);

primaryLowYield = (nShaftIn < max(70, round(0.60 * MaxElecs))) || ...
                  (nLegIn > 0 && nShaftIn < 0.72 * nLegIn);
if primaryLowYield
    support(isLegacy) = support(isLegacy) + 0.10 + 0.10 * legacyChainScore(isLegacy);
    support(weakLegacy) = support(weakLegacy) + 0.12 + 0.10 * legacyChainScore(weakLegacy);
end

baseConf = 0.36 * srcScore + 0.42 * iScore + 0.22 * bScore;
rawConf = baseConf;
rawConf(isLegacy) = baseConf(isLegacy) .* support(isLegacy);

rawKeepThr = 0.82;
if primaryLowYield
    rawKeepThr = 0.74;
end
keepCand = (src == 1) | nearShaft | legacyRescue | (rawConf >= rawKeepThr);
if hasBrainSurf && ~isempty(candMM)
    inCand = LeG_intriangulation(BrainV, BrainF, candMM) == 1;
    if any(inCand)
        D = pdist2(candMM, candMM(inCand,:));
        nearInRad = 7.0;
        strongThr = 0.95;
        if primaryLowYield
            nearInRad = 9.0;
            strongThr = 0.88;
        end
        nearIn = min(D, [], 2) <= nearInRad;
        keepCand = keepCand & (inCand | nearIn | (rawConf >= strongThr) | (src == 1));
    end
end

nLegacyNearShaft = sum(isLegacy & nearShaft);
nLegacyRescue = sum(isLegacy & legacyRescue);
candMM = candMM(keepCand, :);
rawConf = rawConf(keepCand);
src = src(keepCand);
nearShaft = nearShaft(keepCand);
legacyRescue = legacyRescue(keepCand);
fprintf('[AutoElec-Hybrid] Legacy support: near-shaft=%d chain-rescue=%d keptLegacy=%d\n', ...
    nLegacyNearShaft, nLegacyRescue, sum(src == 2));
if isempty(candMM)
    WC = zeros(0,3);
    outConf = zeros(0,1);
    return;
end

if size(candMM,1) > 1
    Z = linkage(candMM, 'average');
    cID = cluster(Z, 'cutoff', 1.8, 'criterion', 'distance');
else
    cID = 1;
end

nG = max(cID);
mm = zeros(nG,3);
cf = zeros(nG,1);
gBoth = false(nG,1);
gLegacyOnly = false(nG,1);
gSupport = zeros(nG,1);
for g = 1:nG
    m = (cID == g);
    w = rawConf(m) + 1e-3;
    w = w / sum(w);
    mm(g,:) = w' * candMM(m,:);
    cf(g) = max(rawConf(m));
    hasShaft = any(src(m) == 1);
    hasLegacy = any(src(m) == 2);
    gBoth(g) = hasShaft && hasLegacy;
    gLegacyOnly(g) = hasLegacy && ~hasShaft;
    gSupport(g) = max(0.6 * double(nearShaft(m)) + 0.4 * double(legacyRescue(m)));
    if gBoth(g)
        cf(g) = cf(g) + 0.15;
    elseif gLegacyOnly(g)
        cf(g) = cf(g) - 0.08 * (1 - gSupport(g));
    end
end
cf = min(cf, 1.2);
cf = max(cf, 0);

if hasBrainSurf && ~isempty(mm)
    in = LeG_intriangulation(BrainV, BrainF, mm) == 1;
    if any(in)
        D = pdist2(mm, mm(in,:));
        nearInRad = 6.0;
        if primaryLowYield
            nearInRad = 8.0;
        end
        nearIn = min(D, [], 2) <= nearInRad;
        keep = in | (nearIn & cf >= 0.35) | ...
               (gBoth & cf >= 0.45) | ...
               (~gLegacyOnly & cf >= (0.72 - 0.08 * primaryLowYield)) | ...
               (gSupport >= (0.78 - 0.10 * primaryLowYield) & cf >= (0.60 - 0.08 * primaryLowYield));
    else
        keep = (cf >= (0.90 - 0.10 * primaryLowYield)) | gBoth | (~gLegacyOnly & cf >= (0.72 - 0.08 * primaryLowYield));
    end
    mm = mm(keep, :);
    cf = cf(keep);
    gLegacyOnly = gLegacyOnly(keep);
    gSupport = gSupport(keep);
end

if ~isempty(mm)
    weakLegacyOnly = gLegacyOnly & (gSupport < (0.55 - 0.10 * primaryLowYield)) & (cf < (0.85 - 0.08 * primaryLowYield));
    keep2 = ~weakLegacyOnly;
    mm = mm(keep2, :);
    cf = cf(keep2);
end

if isempty(mm)
    WC = zeros(0,3);
    outConf = zeros(0,1);
    return;
end

[mm, cf] = pruneWeakOutsideUsingIntensity(mm, cf, Img, MRInfo.mat, p99_8);
[mm, cf] = pruneHybridSparseFragments(mm, cf, BrainV, BrainF, hasBrainSurf);
[mm, cf, ~] = dedupContactsByDistance(mm, cf, ones(size(cf)), 2.0);

[~, ord] = sort(cf, 'descend');
mm = mm(ord, :);
cf = cf(ord);
if size(mm,1) > MaxElecs
    mm = mm(1:MaxElecs, :);
    cf = cf(1:MaxElecs);
end

voxH = MRInfo.mat \ [mm, ones(size(mm,1),1)]';
WC = round(voxH(1:3, :)');
outConf = cf;
end

function [keep, score] = legacyChainSupport(mm, BrainV, BrainF, hasBrainSurf)
n = size(mm,1);
keep = false(n,1);
score = zeros(n,1);
if n < 3
    return;
end

D = pdist2(mm, mm);
A = (D <= 5.5) & (D > 0);
try
    G = graph(A);
    cID = conncomp(G);
catch
    cID = ones(1,n);
end

u = unique(cID);
for i = 1:numel(u)
    idx = find(cID == u(i));
    if numel(idx) < 3
        continue;
    end
    pts = mm(idx, :);
    [dir, mu, span] = shaftModel(pts);
    sv = svd(pts - mu, 'econ');
    if isempty(sv)
        continue;
    end
    if numel(sv) >= 2
        linearity = 1 - sv(2) / max(sv(1), 1e-6);
    else
        linearity = 1;
    end
    linearity = min(max(linearity, 0), 1);

    t = sort((pts - mu) * dir);
    gaps = diff(t);
    gaps = gaps(gaps > 1.0 & gaps < 9.0);
    if isempty(gaps)
        gapMed = 9.0;
        gapReg = 0;
    else
        gapMed = median(gaps);
        gapReg = exp(-abs(gapMed - 3.5) / 1.8) * exp(-std(gaps) / max(1.2, gapMed));
    end

    if hasBrainSurf
        inFrac = mean(LeG_intriangulation(BrainV, BrainF, pts) == 1);
    else
        inFrac = 0.65;
    end

    sizeScore = min(span / 22, 1);
    cntScore = min(log1p(numel(idx)) / log(9), 1);
    compScore = 0.34 * linearity + 0.26 * gapReg + 0.20 * sizeScore + ...
                0.12 * inFrac + 0.08 * cntScore;
    score(idx) = compScore;

    if numel(gaps) >= 2
        spacingOK = (gapMed >= 2.0) && (gapMed <= 6.8) && (std(gaps) <= 2.2);
    elseif numel(gaps) == 1
        spacingOK = (gaps(1) >= 2.0) && (gaps(1) <= 7.0);
    else
        spacingOK = false;
    end

    rescue = (span >= 8.0) && (linearity >= 0.70) && spacingOK && ...
             (~hasBrainSurf || inFrac >= 0.08 || span >= 18.0);
    keep(idx) = keep(idx) | rescue | (compScore >= 0.74 && span >= 10.0);
end
end

function [mm, cf] = pruneHybridSparseFragments(mm, cf, BrainV, BrainF, hasBrainSurf)
if size(mm,1) < 6
    return;
end

D = pdist2(mm, mm);
A = (D <= 6.0) & (D > 0);
try
    G = graph(A);
    cID = conncomp(G);
catch
    cID = ones(1,size(mm,1));
end

keep = true(size(cf));
u = unique(cID);
for i = 1:numel(u)
    idx = find(cID == u(i));
    if numel(idx) >= 4
        continue;
    end
    pts = mm(idx, :);
    [dir, mu, span] = shaftModel(pts);
    if numel(idx) > 1
        t = sort((pts - mu) * dir);
        gaps = diff(t);
        gapMed = median(gaps);
    else
        gapMed = inf;
    end

    if hasBrainSurf
        inFrac = mean(LeG_intriangulation(BrainV, BrainF, pts) == 1);
    else
        inFrac = 0.6;
    end
    weak = (mean(cf(idx)) < 0.72) && (span < 11.0) && (gapMed > 7.0) && (inFrac < 0.25);
    if weak
        keep(idx) = false;
    end
end

mm = mm(keep, :);
cf = cf(keep);
end

function [locI, locB] = sampleLocalAtVox(vox, Img, blob)
imSz = size(Img);
[di,dj,dk] = ndgrid(-1:1, -1:1, -1:1);
offs = [di(:), dj(:), dk(:)];
n = size(vox,1);
locI = zeros(n,1);
locB = zeros(n,1);
for i = 1:n
    rc = round(vox(i,:));
    rc = max(min(rc, imSz), 1);
    nbr = rc + offs;
    nbr = max(nbr, 1);
    nbr(:,1) = min(nbr(:,1), imSz(1));
    nbr(:,2) = min(nbr(:,2), imSz(2));
    nbr(:,3) = min(nbr(:,3), imSz(3));
    idx = sub2ind(imSz, nbr(:,1), nbr(:,2), nbr(:,3));
    locI(i) = mean(Img(idx));
    locB(i) = mean(blob(idx));
end
end

function [mm, cf] = pruneWeakOutsideUsingIntensity(mm, cf, Img, mat, hiThr)
if isempty(mm)
    return;
end
imSz = size(Img);
voxH = mat \ [mm, ones(size(mm,1),1)]';
vox = round(voxH(1:3, :)');
vox = max(vox, 1);
vox(:,1) = min(vox(:,1), imSz(1));
vox(:,2) = min(vox(:,2), imSz(2));
vox(:,3) = min(vox(:,3), imSz(3));
idx = sub2ind(imSz, vox(:,1), vox(:,2), vox(:,3));
v = Img(idx);

keep = true(size(cf));
if numel(cf) > 20
    lowScore = cf < prctile(cf, 25);
    lowInt = v < hiThr;
    keep(lowScore & lowInt) = false;
end
mm = mm(keep, :);
cf = cf(keep);
end


% =========================================================================
%  Shaft-first pipeline (tube/shaft detect first, then contact fit)
% =========================================================================
function [WC, T] = runShaftFirstPipeline(app, Img, MRInfo, ProjSurfRaw, ...
                                         voxDim, ThrR, ThrHU, TMax, TMin, StartTime)
imSz = size(Img);
MaxElecs = getDetectOpt(app, 'maxelecs', 250, [10, 5000]);

vValid = Img(Img > 0);
if isempty(vValid)
    WC = []; T = 0;
    diagPlot(app, ThrHU, zeros(size(ThrHU)), 1, toc(StartTime));
    return;
end

% Robust percentiles for metal/edge modeling.
pct = prctile(vValid(:), [90 95 97 98 99 99.5 99.8 99.9]);
p90 = pct(1); p95 = pct(2);

% Tunables (optional app.detect.*)
corePct          = getDetectOpt(app, 'shaftCorePct', 99.5, [97.5, 99.98]);
haloPct          = getDetectOpt(app, 'shaftHaloPct', 98.6, [95.0, 99.9]);
connRadiusMm     = getDetectOpt(app, 'shaftConnRadiusMm', 1.4, [0.5, 3.0]);
minTubeVolMm3    = getDetectOpt(app, 'shaftMinTubeVolMm3', 1.5, [0.2, 20]);
maxTubeVolMm3    = getDetectOpt(app, 'shaftMaxTubeVolMm3', 420, [20, 3000]);
minShaftLenMm    = getDetectOpt(app, 'shaftMinLenMm', 7.0, [4, 60]);
maxShaftRadiusMm = getDetectOpt(app, 'shaftMaxRadiusMm', 3.2, [1.0, 8.0]);
maxTortuosity    = getDetectOpt(app, 'shaftMaxTortuosity', 1.55, [1.05, 2.5]);
minInsideFrac    = getDetectOpt(app, 'shaftMinInsideFrac', 0.12, [0, 1]);
outsideKeepMm    = getDetectOpt(app, 'shaftOutsideKeepMm', 7.0, [0, 20]);
outsideMinConf   = getDetectOpt(app, 'shaftOutsideMinConf', 0.30, [0, 1]);
spacingDefaultMm = getDetectOpt(app, 'shaftSpacingMm', 3.5, [2.0, 8.0]);
spacingMinMm     = getDetectOpt(app, 'shaftSpacingMinMm', 2.5, [1.5, 6.0]);
spacingMaxMm     = getDetectOpt(app, 'shaftSpacingMaxMm', 6.0, [3.0, 10.0]);
minContactsShaft = round(getDetectOpt(app, 'shaftMinContacts', 2, [1, 12]));
minInterMm       = getDetectOpt(app, 'shaftMinInterContactMm', 2.0, [1.0, 6.0]);

voxVol = prod(voxDim);
minTubeVox = max(3, round(minTubeVolMm3 / max(voxVol, 1e-6)));
maxTubeVox = max(minTubeVox + 1, round(maxTubeVolMm3 / max(voxVol, 1e-6)));

coreThr = prctile(vValid(:), corePct);
haloThr = prctile(vValid(:), min(haloPct, corePct - 0.05));

% If CT is saturated (e.g., both percentiles collapse to 1.0), force separation.
vRange = max(vValid) - min(vValid);
minSep = max(0.01 * vRange, 1e-4);
if (coreThr - haloThr) < minSep
    coreThr = prctile(vValid(:), max(98.8, corePct - 0.4));
    haloThr = prctile(vValid(:), max(95.5, haloPct - 1.2));
end
if coreThr <= haloThr
    coreThr = prctile(vValid(:), 99.0);
    haloThr = prctile(vValid(:), 96.0);
end

fprintf('\n[AutoElec-Shaft] === Shaft-first detection ===\n');
fprintf('[AutoElec-Shaft] core p%.2f=%.1f, halo p%.2f=%.1f, connRadius=%.1f mm\n', ...
    corePct, coreThr, haloPct, haloThr, connRadiusMm);
blob = computeBlobScoreMap(Img, voxDim);
[BrainV, BrainF] = loadBrainSurfaceForAutoElecs(app, ProjSurfRaw);
hasBrainSurf = ~isempty(BrainV) && ~isempty(BrainF);

% Pass 1: stricter.
[finalMM, finalConf, finalShaftID, nShafts] = collectContactsFromTubeMask( ...
    Img, blob, MRInfo, imSz, voxDim, BrainV, BrainF, hasBrainSurf, ...
    coreThr, haloThr, connRadiusMm, ...
    minTubeVox, maxTubeVox, minShaftLenMm, maxTortuosity, maxShaftRadiusMm, minInsideFrac, ...
    outsideKeepMm, outsideMinConf, minContactsShaft, ...
    p90, p95, spacingDefaultMm, spacingMinMm, spacingMaxMm, ...
    0, '[AutoElec-Shaft]');

% Pass 2: automatic relaxed rescue if yield is unexpectedly low.
targetMinContacts = max(40, round(0.35 * MaxElecs));
if size(finalMM,1) < targetMinContacts
    coreThr2 = prctile(vValid(:), max(97.8, corePct - 1.3));
    haloThr2 = prctile(vValid(:), max(95.0, haloPct - 2.0));
    if coreThr2 <= haloThr2
        coreThr2 = prctile(vValid(:), 98.5);
        haloThr2 = prctile(vValid(:), 95.0);
    end
    fprintf('[AutoElec-Shaft] Low yield (%d contacts). Running relaxed pass (core=%.1f halo=%.1f)\n', ...
        size(finalMM,1), coreThr2, haloThr2);

    [mm2, conf2, sid2, nShafts2] = collectContactsFromTubeMask( ...
        Img, blob, MRInfo, imSz, voxDim, BrainV, BrainF, hasBrainSurf, ...
        coreThr2, haloThr2, min(connRadiusMm + 0.5, 2.5), ...
        max(3, round(0.6 * minTubeVox)), round(1.6 * maxTubeVox), ...
        max(5.0, 0.8 * minShaftLenMm), ...
        min(1.9, maxTortuosity + 0.25), ...
        maxShaftRadiusMm + 0.9, ...
        max(0.05, 0.5 * minInsideFrac), ...
        outsideKeepMm + 2.0, max(0.15, outsideMinConf - 0.12), ...
        max(2, minContactsShaft - 1), ...
        p90, p95, spacingDefaultMm, spacingMinMm, spacingMaxMm, ...
        nShafts, '[AutoElec-Shaft-Rlx]');

    if ~isempty(mm2)
        conf2 = 0.85 * conf2;
        finalMM = [finalMM; mm2];
        finalConf = [finalConf; conf2];
        finalShaftID = [finalShaftID; sid2];
        nShafts = nShafts + nShafts2;
    end
end

if ~isempty(finalMM)
    [finalMM, finalConf, finalShaftID] = refineCombinedShafts( ...
        finalMM, finalConf, finalShaftID, ...
        BrainV, BrainF, hasBrainSurf, ...
        spacingMinMm, spacingMaxMm, outsideKeepMm, MaxElecs);
end

if isempty(finalMM)
    WC = []; T = (TMax + TMin) / 2;
    NumObj = zeros(numel(ThrR), 1);
    for kt = 1:numel(ThrR)
        cc0 = bwconncomp(Img > ThrR(kt), 26);
        if cc0.NumObjects == 0, continue; end
        cs0 = cellfun(@numel, cc0.PixelIdxList);
        NumObj(kt) = sum(cs0 >= minTubeVox & cs0 <= maxTubeVox);
    end
    diagPlot(app, ThrHU, NumObj, round(numel(ThrR)/2), toc(StartTime));
    fprintf('[AutoElec-Shaft] No shafts survived geometric/brain filters.\n');
    return;
end

[finalMM, finalConf, finalShaftID] = dedupContactsByDistance(finalMM, finalConf, finalShaftID, minInterMm);
[~, ord] = sort(finalConf, 'descend');
finalMM = finalMM(ord, :);
finalConf = finalConf(ord);
finalShaftID = finalShaftID(ord);

if size(finalMM,1) > MaxElecs
    finalMM = finalMM(1:MaxElecs, :);
    finalConf = finalConf(1:MaxElecs);
    finalShaftID = finalShaftID(1:MaxElecs);
end

voxH = MRInfo.mat \ [finalMM, ones(size(finalMM,1),1)]';
WC = round(voxH(1:3, :)');
T = (TMax + TMin) / 2;

NumObj = zeros(numel(ThrR), 1);
for kt = 1:numel(ThrR)
    cc0 = bwconncomp(Img > ThrR(kt), 26);
    if cc0.NumObjects == 0, continue; end
    cs0 = cellfun(@numel, cc0.PixelIdxList);
    NumObj(kt) = sum(cs0 >= minTubeVox & cs0 <= maxTubeVox);
end
[~, idxD] = min(abs(NumObj - size(WC,1)));
diagPlot(app, ThrHU, NumObj, idxD, toc(StartTime));

fprintf('[AutoElec-Shaft] FINAL: %d shafts, %d contacts in %.1f s\n', ...
    max(finalShaftID), size(WC,1), toc(StartTime));
end


% =========================================================================
%  Helpers for shaft-first pipeline
% =========================================================================
function [finalMM, finalConf, finalShaftID, nShafts] = collectContactsFromTubeMask( ...
    Img, blob, MRInfo, imSz, voxDim, BrainV, BrainF, hasBrainSurf, ...
    coreThr, haloThr, connRadiusMm, ...
    minTubeVox, maxTubeVox, minShaftLenMm, maxTortuosity, maxShaftRadiusMm, minInsideFrac, ...
    outsideKeepMm, outsideMinConf, minContactsShaft, ...
    p90, p95, spacingDefaultMm, spacingMinMm, spacingMaxMm, ...
    shaftIdOffset, logPrefix)

finalMM = zeros(0,3);
finalConf = zeros(0,1);
finalShaftID = zeros(0,1);
nShafts = 0;

maskCore = Img >= coreThr;
maskHalo = Img >= haloThr;
rVox = max(1, round(connRadiusMm / max(min(voxDim), 1e-6)));
se = makeSphereNhood(rVox);
maskTube = imdilate(maskCore, se) & maskHalo;
maskTube = bwareaopen(maskTube, max(3, minTubeVox), 26);

fprintf('%s Pass: core=%.3f halo=%.3f r=%.2f mm coreVox=%d haloVox=%d tubeVox=%d\n', ...
    logPrefix, coreThr, haloThr, connRadiusMm, nnz(maskCore), nnz(maskHalo), nnz(maskTube));

CC = bwconncomp(maskTube, 26);
if CC.NumObjects == 0
    fprintf('%s   no tube components.\n', logPrefix);
    return;
end

nTooSmall = 0;
nTooLarge = 0;
nGeomRejected = 0;
nSplit = 0;

for ci = 1:CC.NumObjects
    idx0 = CC.PixelIdxList{ci};
    nComp0 = numel(idx0);
    if nComp0 < minTubeVox
        nTooSmall = nTooSmall + 1;
        continue;
    end

    compMask0 = false(imSz);
    compMask0(idx0) = true;
    subMasks = {compMask0};
    if nComp0 > maxTubeVox
        % Attempt to split merged tube components by 1-voxel erosion.
        subMasks = splitLargeTubeComponent(compMask0, 40);
        if numel(subMasks) > 1
            nSplit = nSplit + 1;
        elseif nComp0 > (6 * maxTubeVox)
            nTooLarge = nTooLarge + 1;
            continue;
        end
    end

    for sm = 1:numel(subMasks)
        compMask = subMasks{sm};
        idx = find(compMask);
        nComp = numel(idx);
        if nComp < minTubeVox
            continue;
        end

        [~, pathMM, pathS] = extractCenterlineFromComponent(compMask, MRInfo.mat);
        if size(pathMM, 1) < 4
            nGeomRejected = nGeomRejected + 1;
            continue;
        end

        shaftLen = pathS(end) - pathS(1);
        if shaftLen < minShaftLenMm
            nGeomRejected = nGeomRejected + 1;
            continue;
        end

        chordLen = norm(pathMM(end,:) - pathMM(1,:));
        tortuosity = shaftLen / max(chordLen, 1e-6);
        if tortuosity > maxTortuosity
            nGeomRejected = nGeomRejected + 1;
            continue;
        end

        [ii,jj,kk] = ind2sub(imSz, idx);
        compVox = [ii, jj, kk];
        if size(compVox,1) > 1500
            sidx = round(linspace(1, size(compVox,1), 1500));
            compVox = compVox(sidx, :);
        end
        compMM = voxToMM(compVox, MRInfo.mat);
        mu = mean(pathMM, 1);
        [~,~,Vpath] = svd(pathMM - mu, 'econ');
        dir1 = Vpath(:,1);
        tComp = (compMM - mu) * dir1;
        perpDist = sqrt(max(sum((compMM - mu).^2, 2) - tComp.^2, 0));
        radius95 = prctile(perpDist, 95);
        if radius95 > maxShaftRadiusMm
            nGeomRejected = nGeomRejected + 1;
            continue;
        end

        if hasBrainSurf
            sidx = round(linspace(1, size(pathMM,1), min(80, size(pathMM,1))));
            inPath = LeG_intriangulation(BrainV, BrainF, pathMM(sidx,:)) == 1;
            if mean(inPath) < minInsideFrac
                nGeomRejected = nGeomRejected + 1;
                continue;
            end
        end

        [cMM, cConf, cS] = fitContactsAlongShaft(pathMM, pathS, Img, blob, ...
            MRInfo.mat, imSz, p90, p95, spacingDefaultMm, spacingMinMm, spacingMaxMm);
        if isempty(cMM)
            continue;
        end

        if hasBrainSurf
            inC = LeG_intriangulation(BrainV, BrainF, cMM) == 1;
            if any(inC)
                sInside = cS(inC);
                sMin = min(sInside) - outsideKeepMm;
                sMax = max(sInside) + outsideKeepMm;
                keepC = inC | ((cS >= sMin & cS <= sMax) & cConf >= outsideMinConf);
                cMM = cMM(keepC, :);
                cConf = cConf(keepC);
                cS = cS(keepC);
            else
                cMM = zeros(0,3);
                cConf = zeros(0,1);
            end
        end

        if size(cMM,1) < minContactsShaft
            continue;
        end

        nShafts = nShafts + 1;
        sid = shaftIdOffset + nShafts;
        finalMM = [finalMM; cMM]; %#ok<AGROW>
        finalConf = [finalConf; cConf(:)]; %#ok<AGROW>
        finalShaftID = [finalShaftID; sid * ones(size(cMM,1),1)]; %#ok<AGROW>

        fprintf('%s   Shaft %d: L=%.1f mm, tort=%.2f, rad95=%.2f mm -> %d contacts\n', ...
            logPrefix, sid, shaftLen, tortuosity, radius95, size(cMM,1));
    end
end

fprintf('%s Summary: CC=%d, splitCC=%d, tooSmall=%d, tooLarge=%d, geomReject=%d, shafts=%d, contacts=%d\n', ...
    logPrefix, CC.NumObjects, nSplit, nTooSmall, nTooLarge, nGeomRejected, nShafts, size(finalMM,1));
end

function val = getDetectOpt(app, fieldName, defaultVal, bounds)
val = defaultVal;
det = getDetectStruct(app);
if isstruct(det) && isfield(det, fieldName) && ~isempty(det.(fieldName))
    raw = det.(fieldName);
    if isnumeric(raw) && isfinite(raw(1))
        val = double(raw(1));
    end
end
if nargin >= 4 && numel(bounds) == 2
    val = min(max(val, bounds(1)), bounds(2));
end
end

function det = getDetectStruct(app)
det = struct;
try
    if isstruct(app)
        if isfield(app, 'detect') && isstruct(app.detect)
            det = app.detect;
        end
    elseif isobject(app)
        if isprop(app, 'detect') && isstruct(app.detect)
            det = app.detect;
        end
    end
catch
    det = struct;
end
end

function se = makeSphereNhood(r)
[ii,jj,kk] = ndgrid(-r:r, -r:r, -r:r);
se = (ii.^2 + jj.^2 + kk.^2) <= (r^2);
end

function blob = computeBlobScoreMap(Img, voxDim)
imSz = size(Img);
sigmas_mm = [0.4 0.8 1.2];
blob = zeros(imSz, 'single');
for s = 1:numel(sigmas_mm)
    sigma_mm = sigmas_mm(s);
    sigma_vox = sigma_mm ./ voxDim(:)';
    g1 = imgaussfilt3(Img, sigma_vox);
    g2 = imgaussfilt3(Img, sigma_vox * 1.6);
    blob = max(blob, -(g2 - g1) * sigma_mm^2);
end
bClip = prctile(blob(blob > 0), 99.5);
if ~isempty(bClip) && bClip > 0
    blob = blob / bClip;
end
end

function [BrainV, BrainF] = loadBrainSurfaceForAutoElecs(app, ProjSurfRaw)
BrainV = [];
BrainF = [];
if isfield(app, 'SurfacesFile') && ~isempty(app.SurfacesFile)
    try
        loaded = load(app.SurfacesFile);
        if isfield(loaded, 'BrainSurfRaw')
            S = loaded.BrainSurfRaw;
        elseif isfield(loaded, 'k') && isfield(loaded.k, 'BrainSurfRaw')
            S = loaded.k.BrainSurfRaw;
        else
            fns = fieldnames(loaded);
            if ~isempty(fns) && isfield(loaded.(fns{1}), 'BrainSurfRaw')
                S = loaded.(fns{1}).BrainSurfRaw;
            else
                S = [];
            end
        end
        if ~isempty(S) && isfield(S, 'vertices') && isfield(S, 'faces')
            BrainV = S.vertices;
            BrainF = S.faces;
            return;
        end
    catch
    end
end
if ~isempty(ProjSurfRaw) && isfield(ProjSurfRaw, 'vertices') && isfield(ProjSurfRaw, 'faces')
    BrainV = ProjSurfRaw.vertices;
    BrainF = ProjSurfRaw.faces;
end
end

function mm = voxToMM(vox, mat)
mmH = [vox, ones(size(vox,1),1)] * mat';
mm = mmH(:,1:3);
end

function [pathVox, pathMM, sPath] = extractCenterlineFromComponent(compMask, mat)
pathVox = zeros(0,3);
pathMM = zeros(0,3);
sPath = zeros(0,1);

try
    sk = bwskel(compMask, 'MinBranchLength', 2);
catch
    try
        sk = bwmorph3(compMask, 'skel', Inf);
    catch
        sk = compMask;
    end
end
if ~any(sk(:))
    sk = compMask;
end

idx = find(sk);
if numel(idx) < 2
    return;
end
[ii,jj,kk] = ind2sub(size(sk), idx);
pts = [ii, jj, kk];

pathVox = orderPointsByGraphPath(pts);
if size(pathVox,1) < 2
    return;
end

pathMM = voxToMM(pathVox, mat);
if size(pathMM,1) > 4
    pathMM(2:end-1,:) = (pathMM(1:end-2,:) + pathMM(2:end-1,:) + pathMM(3:end,:)) / 3;
end

d = sqrt(sum(diff(pathMM,1,1).^2,2));
keep = [true; d > 1e-3];
pathMM = pathMM(keep, :);
if size(pathMM,1) < 2
    pathVox = zeros(0,3);
    sPath = zeros(0,1);
    return;
end
sPath = [0; cumsum(sqrt(sum(diff(pathMM,1,1).^2,2)))];
end

function path = orderPointsByGraphPath(pts)
path = zeros(0,3);
n = size(pts,1);
if n < 2
    path = pts;
    return;
end

if n > 2500
    mu = mean(pts,1);
    [~,~,V] = svd(pts - mu, 'econ');
    t = (pts - mu) * V(:,1);
    [~, ord] = sort(t);
    path = pts(ord, :);
    return;
end

try
    D = pdist2(pts, pts, 'chebychev');
    A = (D <= 1) & (D > 0);
    G = graph(A);

    compId = conncomp(G);
    if numel(unique(compId)) > 1
        compCounts = accumarray(compId(:), 1);
        [~, mx] = max(compCounts);
        keep = (compId == mx);
        pts = pts(keep, :);
        D = pdist2(pts, pts, 'chebychev');
        A = (D <= 1) & (D > 0);
        G = graph(A);
    end

    deg = degree(G);
    endPts = find(deg == 1);
    if numel(endPts) >= 2
        DE = distances(G, endPts, endPts);
        DE(~isfinite(DE)) = -Inf;
        [mxVal, mxIdx] = max(DE(:));
        if isfinite(mxVal) && mxVal > 0
            [e1, e2] = ind2sub(size(DE), mxIdx);
            pIdx = shortestpath(G, endPts(e1), endPts(e2));
            path = pts(pIdx, :);
            return;
        end
    end
catch
end

mu = mean(pts,1);
[~,~,V] = svd(pts - mu, 'econ');
t = (pts - mu) * V(:,1);
[~, ord] = sort(t);
path = pts(ord, :);
end

function [contactMM, contactConf, contactS] = fitContactsAlongShaft(pathMM, pathS, Img, blob, ...
    mat, imSz, p90, p95, spacingDefaultMm, spacingMinMm, spacingMaxMm)

contactMM = zeros(0,3);
contactConf = zeros(0,1);
contactS = zeros(0,1);

if size(pathMM,1) < 3
    return;
end

ds = 0.5;
[sampMM, sampS] = resamplePolyline(pathMM, ds);
if size(sampMM,1) < 5
    return;
end

[localI, localB] = sampleProfileOnPath(sampMM, Img, blob, mat, imSz);
iNorm = max(0, (localI - p90) / max(p95 - p90, 1e-6));
bNorm = max(0, localB / max(prctile(localB, 95), 1e-6));
score = 0.7 * iNorm + 0.3 * bNorm;
if max(score) < 1e-6
    score = localI / max(max(localI), 1e-6);
end

minPeakDist = max(1, round(2.2 / ds));
peakThr = max(prctile(score, 55), 0.05);
peakIdx = findLocalMaxima1D(score, minPeakDist, peakThr);

if numel(peakIdx) >= 2
    spacingMm = estimateSpacing(sampS(peakIdx), spacingMinMm, spacingMaxMm, spacingDefaultMm);
else
    spacingMm = spacingDefaultMm;
end

[gridS, gridScore, peakDist] = fitGridContacts(sampS, score, peakIdx, spacingMm, ds);
if isempty(gridS)
    return;
end

keep = gridScore >= max(prctile(score, 45), 0.04) | peakDist <= 2.0;
contactS = gridS(keep);
if isempty(contactS)
    contactS = sampS(peakIdx);
end
% Ensure minimum density implied by shaft length and allowed max spacing.
minExpected = max(2, floor((sampS(end) - sampS(1)) / max(spacingMaxMm, 1e-3)) + 1);
if numel(contactS) < minExpected
    contactS = [contactS(:); gridS(:)]; %#ok<AGROW>
end
contactS = mergeClosePositions(contactS(:), 1.0);
if isempty(contactS)
    return;
end

contactMM = [interp1(sampS, sampMM(:,1), contactS, 'linear', 'extrap'), ...
             interp1(sampS, sampMM(:,2), contactS, 'linear', 'extrap'), ...
             interp1(sampS, sampMM(:,3), contactS, 'linear', 'extrap')];

scoreAt = interp1(sampS, score, contactS, 'linear', 0);
if ~isempty(peakIdx)
    peakS = sampS(peakIdx(:));
    dMat = abs(bsxfun(@minus, contactS(:), peakS(:)'));
    peakDist = min(dMat, [], 2);
else
    peakDist = inf(size(contactS));
end
contactConf = min(1, max(0.05, 0.65 * scoreAt + 0.35 * exp(-peakDist/1.5)));

maxN = max(2, ceil((pathS(end) - pathS(1)) / max(spacingMinMm, 1.5)) + 2);
if numel(contactConf) > maxN
    [~, o] = sort(contactConf, 'descend');
    keepIdx = o(1:maxN);
    contactMM = contactMM(keepIdx, :);
    contactConf = contactConf(keepIdx);
    contactS = contactS(keepIdx);
end

[contactS, o2] = sort(contactS);
contactMM = contactMM(o2, :);
contactConf = contactConf(o2);
end

function [sampMM, sampS] = resamplePolyline(pathMM, ds)
d = sqrt(sum(diff(pathMM,1,1).^2,2));
s = [0; cumsum(d)];
[s, ia] = unique(s, 'stable');
pathMM = pathMM(ia, :);
if numel(s) < 2
    sampMM = pathMM;
    sampS = s;
    return;
end

sampS = (0:ds:s(end))';
if sampS(end) < s(end)
    sampS = [sampS; s(end)]; %#ok<AGROW>
end
sampMM = [interp1(s, pathMM(:,1), sampS, 'linear'), ...
          interp1(s, pathMM(:,2), sampS, 'linear'), ...
          interp1(s, pathMM(:,3), sampS, 'linear')];
end

function [localI, localB] = sampleProfileOnPath(pathMM, Img, blob, mat, imSz)
n = size(pathMM,1);
localI = zeros(n,1);
localB = zeros(n,1);
[di,dj,dk] = ndgrid(-1:1, -1:1, -1:1);
offs = [di(:), dj(:), dk(:)];
for i = 1:n
    v = mat \ [pathMM(i,:)'; 1];
    rc = round(v(1:3))';
    rc = max(min(rc, imSz), 1);
    nbr = rc + offs;
    nbr = max(nbr, 1);
    nbr(:,1) = min(nbr(:,1), imSz(1));
    nbr(:,2) = min(nbr(:,2), imSz(2));
    nbr(:,3) = min(nbr(:,3), imSz(3));
    idx = sub2ind(imSz, nbr(:,1), nbr(:,2), nbr(:,3));
    localI(i) = mean(Img(idx));
    localB(i) = mean(blob(idx));
end
end

function peakIdx = findLocalMaxima1D(x, minSep, minThr)
peakIdx = find(x(2:end-1) >= x(1:end-2) & x(2:end-1) >= x(3:end)) + 1;
if isempty(peakIdx)
    return;
end
peakIdx = peakIdx(x(peakIdx) >= minThr);
if isempty(peakIdx)
    return;
end

[~, ord] = sort(x(peakIdx), 'descend');
cand = peakIdx(ord);
keep = false(size(x));
sel = zeros(0,1);
for i = 1:numel(cand)
    c = cand(i);
    lo = max(1, c - minSep);
    hi = min(numel(x), c + minSep);
    if ~any(keep(lo:hi))
        sel(end+1,1) = c; %#ok<AGROW>
        keep(c) = true;
    end
end
peakIdx = sort(sel);
end

function spacingMm = estimateSpacing(peakS, spacingMinMm, spacingMaxMm, defaultSpacing)
spacingMm = defaultSpacing;
if numel(peakS) < 2
    return;
end
gaps = diff(peakS(:));
gaps = gaps(gaps > 1.0);
if isempty(gaps)
    return;
end

cand = spacingMinMm:0.1:spacingMaxMm;
scores = zeros(size(cand));
for i = 1:numel(cand)
    r = gaps / cand(i);
    k = max(round(r), 1);
    err = abs(r - k);
    scores(i) = sum(exp(-(err/0.22).^2)) - 0.1*sum(err);
end
[~, bi] = max(scores);
spacingMm = cand(bi);
end

function [gridS, gridScore, peakDist] = fitGridContacts(sampS, score, peakIdx, spacingMm, ds)
gridS = zeros(0,1);
gridScore = zeros(0,1);
peakDist = zeros(0,1);
if sampS(end) <= 0 || spacingMm <= 0
    return;
end

offCand = (0:ds:spacingMm)';
bestVal = -Inf;
bestGrid = [];
bestScore = [];
bestDist = [];

if isempty(peakIdx)
    peakS = [];
else
    peakS = sampS(peakIdx(:));
end

for oi = 1:numel(offCand)
    g = (offCand(oi):spacingMm:sampS(end))';
    if numel(g) < 2, continue; end
    vals = interp1(sampS, score, g, 'linear', 0);
    if isempty(peakS)
        dMin = inf(size(g));
        align = 0;
    else
        dMat = abs(bsxfun(@minus, g(:), peakS(:)'));
        dMin = min(dMat, [], 2);
        align = mean(exp(-(dMin/1.2).^2));
    end
    val = mean(vals) + 0.35 * align;
    if val > bestVal
        bestVal = val;
        bestGrid = g;
        bestScore = vals;
        bestDist = dMin;
    end
end

gridS = bestGrid;
gridScore = bestScore;
peakDist = bestDist;
end

function x = mergeClosePositions(x, tol)
if isempty(x)
    return;
end
x = sort(x(:));
y = zeros(0,1);
i = 1;
while i <= numel(x)
    j = i;
    while j < numel(x) && abs(x(j+1) - x(j)) <= tol
        j = j + 1;
    end
    y(end+1,1) = mean(x(i:j)); %#ok<AGROW>
    i = j + 1;
end
x = y;
end

function subMasks = splitLargeTubeComponent(compMask, maxParts)
subMasks = {compMask};
try
    se1 = makeSphereNhood(1);
    er = imerode(compMask, se1);
    if nnz(er) < 5
        return;
    end
    CCe = bwconncomp(er, 26);
    if CCe.NumObjects <= 1 || CCe.NumObjects > maxParts
        return;
    end
    subMasks = cell(CCe.NumObjects, 1);
    for i = 1:CCe.NumObjects
        m = false(size(compMask));
        m(CCe.PixelIdxList{i}) = true;
        m = imdilate(m, se1) & compMask;
        m = bwareaopen(m, 3, 26);
        subMasks{i} = m;
    end
catch
    subMasks = {compMask};
end
end

function [mmOut, confOut, sidOut] = refineCombinedShafts( ...
    mm, conf, sid, BrainV, BrainF, hasBrainSurf, spacingMinMm, spacingMaxMm, outsideKeepMm, MaxElecs)

if isempty(mm)
    mmOut = mm; confOut = conf; sidOut = sid;
    return;
end

% 1) Merge fragmented shafts that are collinear and spatially overlapping.
[mm, conf, sid] = mergeCollinearShafts(mm, conf, sid);

uS = unique(sid(:))';
nS = numel(uS);
if nS == 0
    mmOut = zeros(0,3); confOut = zeros(0,1); sidOut = zeros(0,1);
    return;
end

scoreS = zeros(nS,1);
nContS = zeros(nS,1);
keepContact = true(size(conf));

for i = 1:nS
    s = uS(i);
    idx = (sid == s);
    pts = mm(idx, :);
    cf  = conf(idx);
    nC = size(pts,1);
    nContS(i) = nC;

    [dir, mu, span] = shaftModel(pts);
    if nC >= 3
        spacingReg = spacingRegularity(pts, dir, mu, spacingMinMm, spacingMaxMm);
    elseif nC == 2
        spacingReg = 0.45;
    else
        spacingReg = 0.2;
    end

    if hasBrainSurf
        in = LeG_intriangulation(BrainV, BrainF, pts) == 1;
        inFrac = mean(in);
        t = (pts - mu) * dir;
        if any(in)
            tIn = t(in);
            tMin = min(tIn) - outsideKeepMm;
            tMax = max(tIn) + outsideKeepMm;
            keepLocal = in | (((t >= tMin) & (t <= tMax)) & (cf >= 0.25));
            ex = max(0, t - tMax) + max(0, tMin - t);
            outsidePenalty = exp(-mean(ex)/2.0);
        else
            keepLocal = (cf >= 0.55) & (span >= 15);
            outsidePenalty = 0.5;
        end
    else
        inFrac = 0.7;
        keepLocal = true(size(cf));
        outsidePenalty = 1.0;
    end

    tmp = find(idx);
    keepContact(tmp(~keepLocal)) = false;
    ptsK = pts(keepLocal, :);
    cfK = cf(keepLocal);
    if size(ptsK,1) < 1
        scoreS(i) = 0;
        nContS(i) = 0;
        continue;
    end

    [~, ~, spanK] = shaftModel(ptsK);
    meanCf = mean(cfK);
    sLen = min(spanK / 30, 1);
    sConf = min(meanCf / 0.55, 1);
    sCount = min(log1p(size(ptsK,1)) / log(8), 1);
    sBrain = inFrac;
    scoreS(i) = outsidePenalty * (0.23*sLen + 0.24*sConf + 0.22*spacingReg + 0.21*sBrain + 0.10*sCount);
end

mm = mm(keepContact, :);
conf = conf(keepContact);
sid = sid(keepContact);

uS = unique(sid(:))';
nS = numel(uS);
if nS == 0
    mmOut = zeros(0,3); confOut = zeros(0,1); sidOut = zeros(0,1);
    return;
end

scoreNow = zeros(nS,1);
nContNow = zeros(nS,1);
for i = 1:nS
    s = uS(i);
    idx = (sid == s);
    pts = mm(idx, :);
    cf = conf(idx);
    nContNow(i) = size(pts,1);
    [dir, mu, span] = shaftModel(pts);

    if nContNow(i) >= 3
        spacingReg = spacingRegularity(pts, dir, mu, spacingMinMm, spacingMaxMm);
    elseif nContNow(i) == 2
        spacingReg = 0.45;
    else
        spacingReg = 0.2;
    end

    if hasBrainSurf
        in = LeG_intriangulation(BrainV, BrainF, pts) == 1;
        inFrac = mean(in);
        t = (pts - mu) * dir;
        if any(in)
            tIn = t(in);
            tMin = min(tIn) - outsideKeepMm;
            tMax = max(tIn) + outsideKeepMm;
            ex = max(0, t - tMax) + max(0, tMin - t);
            outsidePenalty = exp(-mean(ex)/2.0);
        else
            outsidePenalty = 0.5;
        end
    else
        inFrac = 0.7;
        outsidePenalty = 1.0;
    end

    sLen = min(span / 30, 1);
    sConf = min(mean(cf) / 0.55, 1);
    sCount = min(log1p(nContNow(i)) / log(8), 1);
    scoreNow(i) = outsidePenalty * (0.23*sLen + 0.24*sConf + 0.22*spacingReg + 0.21*inFrac + 0.10*sCount);
end

thr = max(0.30, prctile(scoreNow, 30) - 0.05);
keepS = scoreNow >= thr;

% Keep strong short shafts, reject weak short fragments.
for i = 1:nS
    if ~keepS(i), continue; end
    s = uS(i);
    idx = (sid == s);
    pts = mm(idx,:);
    [~, ~, span] = shaftModel(pts);
    if nContNow(i) <= 2 && scoreNow(i) < 0.62 && span < 12
        keepS(i) = false;
    elseif nContNow(i) <= 3 && scoreNow(i) < 0.45 && span < 9
        keepS(i) = false;
    end
end

targetMin = max(70, round(0.45 * MaxElecs));
nKeepCont = sum(ismember(sid, uS(keepS)));
if nKeepCont < targetMin
    thr2 = max(0.20, thr - 0.10);
    keepS = scoreNow >= thr2;
end

keepMask = ismember(sid, uS(keepS));
mmOut = mm(keepMask, :);
confOut = conf(keepMask);
sidOut = sid(keepMask);

fprintf('[AutoElec-Shaft] Refinement: shafts %d -> %d, contacts %d -> %d, scoreThr=%.2f\n', ...
    nS, sum(keepS), size(mm,1), size(mmOut,1), thr);
end

function [mm, conf, sid] = mergeCollinearShafts(mm, conf, sid)
if isempty(mm)
    return;
end
changed = true;
while changed
    changed = false;
    u = unique(sid(:))';
    if numel(u) < 2
        break;
    end
    for i = 1:numel(u)-1
        idx1 = (sid == u(i));
        pts1 = mm(idx1, :);
        if size(pts1,1) < 2, continue; end
        [d1, ~, ~] = shaftModel(pts1);
        for j = i+1:numel(u)
            idx2 = (sid == u(j));
            pts2 = mm(idx2, :);
            if size(pts2,1) < 2, continue; end
            [d2, ~, ~] = shaftModel(pts2);

            ang = acosd(min(1, max(-1, abs(dot(d1, d2)))));
            if ang > 18
                continue;
            end
            ld = lineDistance(pts1, d1, pts2, d2);
            if ld > 2.5
                continue;
            end
            gap = shaftIntervalGap(pts1, pts2, d1, d2);
            if gap > 6
                continue;
            end
            sid(idx2) = u(i);
            changed = true;
            break;
        end
        if changed
            break;
        end
    end
end
end

function [dir, mu, span] = shaftModel(pts)
mu = mean(pts, 1);
if size(pts,1) >= 2
    [~,~,V] = svd(pts - mu, 'econ');
    dir = V(:,1);
else
    dir = [1;0;0];
end
t = (pts - mu) * dir;
span = max(t) - min(t);
end

function s = spacingRegularity(pts, dir, mu, spacingMinMm, spacingMaxMm)
t = sort((pts - mu) * dir);
g = diff(t);
g = g(g > 0.7);
if isempty(g)
    s = 0.3;
    return;
end
cand = spacingMinMm:0.1:spacingMaxMm;
sc = zeros(numel(cand),1);
for i = 1:numel(cand)
    r = g / cand(i);
    k = max(round(r), 1);
    e = abs(r - k);
    sc(i) = mean(exp(-(e/0.25).^2));
end
s = max(sc);
end

function d = lineDistance(pts1, d1, pts2, d2)
mu1 = mean(pts1, 1);
mu2 = mean(pts2, 1);
cr = cross(d1(:), d2(:));
if norm(cr) < 1e-6
    d = norm(cross((mu2 - mu1)', d1(:)));
else
    d = abs(dot((mu2 - mu1), cr')) / norm(cr);
end
end

function gap = shaftIntervalGap(pts1, pts2, d1, d2)
sgn = sign(dot(d1(:), d2(:)));
if sgn == 0, sgn = 1; end
d = d1(:) + sgn*d2(:);
d = d / max(norm(d), 1e-6);
t1 = pts1 * d;
t2 = pts2 * d;
gap = max(0, max(min(t1)-max(t2), min(t2)-max(t1)));
end

function [mm, conf, sid] = dedupContactsByDistance(mm, conf, sid, minDistMm)
if size(mm,1) < 2
    return;
end
[~, ord] = sort(conf, 'descend');
mm = mm(ord, :);
conf = conf(ord);
sid = sid(ord);

keep = true(size(conf));
for i = 1:size(mm,1)
    if ~keep(i), continue; end
    d = sqrt(sum((mm(i+1:end,:) - mm(i,:)).^2, 2));
    rm = find(d < minDistMm) + i;
    keep(rm) = false;
end
mm = mm(keep, :);
conf = conf(keep);
sid = sid(keep);
end


% =========================================================================
%  diagPlot — backward-compatible diagnostic figure
% =========================================================================
function diagPlot(app, ThrHU, NumObj, idx, elapsedTime)
    fH = figure('Position', [50 50 400 400], ...
                'Name',     app.PatientIDStr);
    aH = axes('Parent', fH);
    plot(aH, ThrHU, NumObj, '-');
    hold(aH, 'on');

    if idx >= 1 && idx <= numel(ThrHU)
        plot(aH, ThrHU(idx), NumObj(idx), 'or');
    end

    xlabel(aH, 'Threshold (HU)');
    ylabel(aH, '# Detections');

    TMaxHU = (app.CTRng(4) - app.CTInfo.pinfo(2)) ./ app.CTInfo.pinfo(1);
    TMinHU = (app.CTRng(3) - app.CTInfo.pinfo(2)) ./ app.CTInfo.pinfo(1);
    title(aH, sprintf('%0.1f s  (%0.0f, %0.0f) HU', ...
          elapsedTime, TMinHU, TMaxHU));
end


% =========================================================================
%  diagPlotsAutoElecs — diagnostic figures for sensitivity / scan comparison
% =========================================================================
%  Opens figures: confidence distribution, 3D contacts (colored by confidence/shaft),
%  optional surface overlay, inter-contact distance distribution.
function diagPlotsAutoElecs(app, ThrHU, NumObj, mrgMM, confidence, assigned, ...
    finalMM_afterPhase3, finalConf_afterPhase3, finalShaftID_afterPhase3, ...
    nAfterBoneFilter, finalMM, finalConf, finalShaftID, elapsedTime)
    nC = size(finalMM, 1);
    nP1 = size(mrgMM, 1);
    sz = max(8, 80 - nC);

    % Load brain surface: SurfacesFile first (BrainSurfRaw), then ProjSurfRaw
    surfV = []; surfF = [];
    if isfield(app, 'SurfacesFile') && ~isempty(app.SurfacesFile)
        try
            loaded = load(app.SurfacesFile);
            if isfield(loaded, 'BrainSurfRaw'), S = loaded.BrainSurfRaw;
            elseif isfield(loaded, 'k') && isfield(loaded.k, 'BrainSurfRaw'), S = loaded.k.BrainSurfRaw;
            else, fns = fieldnames(loaded); if ~isempty(fns) && isfield(loaded.(fns{1}), 'BrainSurfRaw'), S = loaded.(fns{1}).BrainSurfRaw; else, S = []; end; end
            if ~isempty(S) && isfield(S, 'vertices') && isfield(S, 'faces')
                surfV = S.vertices; surfF = S.faces;
            end
        catch
        end
    end
    if isempty(surfV) && isfield(app, 'ProjSurfRaw') && ~isempty(app.ProjSurfRaw) && isfield(app.ProjSurfRaw, 'vertices') && isfield(app.ProjSurfRaw, 'faces')
        surfV = app.ProjSurfRaw.vertices;
        surfF = app.ProjSurfRaw.faces;
    end

    % ----- 1) Pipeline funnel: contact count at each stage -----
    figure('Position', [50 50 520 380], 'Name', [app.PatientIDStr ' — Pipeline funnel']);
    nAssigned = sum(assigned);
    nAfterP3 = size(finalMM_afterPhase3, 1);
    counts = [nP1; nAssigned; nAfterP3; nAfterBoneFilter; nC];
    labels = {'Phase 1 (merged)', 'Assigned to shaft', 'After Phase 3', 'After filter', 'Final'};
    bar(counts, 'FaceColor', [0.35 0.55 0.8], 'EdgeColor', 'none');
    set(gca, 'XTickLabel', labels);
    ylabel('Number of contacts');
    title(sprintf('Pipeline: %d -> %d -> %d -> %d -> %d', counts(1), counts(2), counts(3), counts(4), counts(5)));
    grid on;

    % ----- 2) 3D Phase 1: all candidates + surface, colored by confidence -----
    sz1 = max(6, 60 - round(nP1/20));
    figure('Position', [60 60 600 480], 'Name', [app.PatientIDStr ' — Phase 1 candidates']);
    ax = gca;
    if ~isempty(surfV), patch(ax, 'Vertices', surfV, 'Faces', surfF, 'FaceColor', [0.85 0.85 0.92], 'FaceAlpha', 0.4, 'EdgeColor', 'none'); hold(ax, 'on'); end
    scatter3(ax, mrgMM(:,1), mrgMM(:,2), mrgMM(:,3), sz1, confidence, 'filled');
    colormap(ax, 'parula');
    cb = colorbar(ax); cb.Label.String = 'Confidence';
    xlabel(ax, 'X (mm)'); ylabel(ax, 'Y (mm)'); zlabel(ax, 'Z (mm)');
    title(ax, sprintf('Phase 1 merged (n=%d) — colored by confidence', nP1));
    axis(ax, 'equal'); grid(ax, 'on'); view(ax, 3);

    % ----- 3) Phase 1 outcome: shaft / singleton kept / dropped + surface -----
    status = zeros(nP1, 1);
    status(assigned) = 1;
    for i = find(~assigned)'
        if min(pdist2(mrgMM(i,:), finalMM)) <= 3.0, status(i) = 2; else, status(i) = 3; end
    end
    figure('Position', [70 70 600 480], 'Name', [app.PatientIDStr ' — Phase 1 outcome']);
    ax = gca;
    if ~isempty(surfV), patch(ax, 'Vertices', surfV, 'Faces', surfF, 'FaceColor', [0.85 0.85 0.92], 'FaceAlpha', 0.4, 'EdgeColor', 'none'); hold(ax, 'on'); end
    col = [0 0.7 0; 0 0.4 0.9; 0.9 0.25 0.2];
    for k = 1:3
        idx = (status == k);
        if any(idx), scatter3(ax, mrgMM(idx,1), mrgMM(idx,2), mrgMM(idx,3), sz1, col(k,:), 'filled'); end
    end
    legend(ax, 'Shaft', 'Singleton kept', 'Dropped', 'Location', 'best');
    xlabel(ax, 'X (mm)'); ylabel(ax, 'Y (mm)'); zlabel(ax, 'Z (mm)');
    title(ax, sprintf('Phase 1 outcome: %d shaft, %d kept, %d dropped', sum(status==1), sum(status==2), sum(status==3)));
    axis(ax, 'equal'); grid(ax, 'on'); view(ax, 3);

    % ----- 4) 3D After Phase 3 + surface -----
    if nAfterP3 > 0
        figure('Position', [80 80 600 480], 'Name', [app.PatientIDStr ' — After Phase 3']);
        ax = gca;
        if ~isempty(surfV), patch(ax, 'Vertices', surfV, 'Faces', surfF, 'FaceColor', [0.85 0.85 0.92], 'FaceAlpha', 0.4, 'EdgeColor', 'none'); hold(ax, 'on'); end
        scatter3(ax, finalMM_afterPhase3(:,1), finalMM_afterPhase3(:,2), finalMM_afterPhase3(:,3), max(6, sz), finalConf_afterPhase3, 'filled');
        colormap(ax, 'parula'); colorbar(ax);
        xlabel(ax, 'X (mm)'); ylabel(ax, 'Y (mm)'); zlabel(ax, 'Z (mm)');
        title(ax, sprintf('After Phase 3 (n=%d) — before bone/linearity filter', nAfterP3));
        axis(ax, 'equal'); grid(ax, 'on'); view(ax, 3);
    end

    % ----- 4b) Kept vs dropped after Phase 3 (two colors, same figure) -----
    if nAfterP3 > 0 && nC > 0
        matchMm = 5.0;
        distToFinal = pdist2(finalMM_afterPhase3, finalMM);
        kept = min(distToFinal, [], 2) <= matchMm;
        nKept = sum(kept);
        nDropped = sum(~kept);
        figure('Position', [82 82 600 480], 'Name', [app.PatientIDStr ' — Kept vs dropped after Phase 3']);
        ax = gca;
        hold(ax, 'on');
        if ~isempty(surfV), patch(ax, 'Vertices', surfV, 'Faces', surfF, 'FaceColor', [0.85 0.85 0.92], 'FaceAlpha', 0.4, 'EdgeColor', 'none'); end
        szP3 = max(6, 60 - round(nAfterP3/20));
        % Plot both: red = dropped, green = kept. Build legend only for series we add.
        legEntries = {};
        if nDropped > 0
            scatter3(ax, finalMM_afterPhase3(~kept,1), finalMM_afterPhase3(~kept,2), finalMM_afterPhase3(~kept,3), szP3, [0.85 0.2 0.15], 'filled', 'DisplayName', sprintf('Dropped (n=%d)', nDropped));
            legEntries{end+1} = sprintf('Dropped (n=%d)', nDropped);
        end
        if nKept > 0
            scatter3(ax, finalMM_afterPhase3(kept,1), finalMM_afterPhase3(kept,2), finalMM_afterPhase3(kept,3), szP3, [0 0.65 0.2], 'filled', 'DisplayName', sprintf('Kept (n=%d)', nKept));
            legEntries{end+1} = sprintf('Kept (n=%d)', nKept);
        end
        if ~isempty(legEntries), legend(ax, legEntries, 'Location', 'best'); end
        xlabel(ax, 'X (mm)'); ylabel(ax, 'Y (mm)'); zlabel(ax, 'Z (mm)');
        title(ax, 'After Phase 3: kept (green) vs dropped by bone/linearity/Phase 4 (red)');
        axis(ax, 'equal'); grid(ax, 'on'); view(ax, 3);
    end

    % ----- 5) Confidence: Phase 1 vs Final -----
    figure('Position', [85 85 560 420], 'Name', [app.PatientIDStr ' — Confidence (Phase 1 vs Final)']);
    subplot(2,1,1);
    histogram(confidence, 25, 'FaceColor', [0.5 0.6 0.9], 'EdgeColor', 'none');
    xlabel('Confidence'); ylabel('Count');
    title(sprintf('Phase 1 merged (n=%d) — mean %.2f', nP1, mean(confidence)));
    grid on;
    subplot(2,1,2);
    histogram(finalConf, 25, 'FaceColor', [0.2 0.5 0.8], 'EdgeColor', 'none');
    xlabel('Confidence'); ylabel('Count');
    title(sprintf('Final (n=%d) — mean %.2f, median %.2f', nC, mean(finalConf), median(finalConf)));
    grid on;

    % ----- 6) 3D Final + surface, colored by confidence -----
    figure('Position', [95 95 600 480], 'Name', [app.PatientIDStr ' — Final contacts + surface']);
    ax = gca;
    if ~isempty(surfV), patch(ax, 'Vertices', surfV, 'Faces', surfF, 'FaceColor', [0.85 0.85 0.92], 'FaceAlpha', 0.4, 'EdgeColor', 'none'); hold(ax, 'on'); end
    scatter3(ax, finalMM(:,1), finalMM(:,2), finalMM(:,3), sz, finalConf, 'filled');
    colormap(ax, 'parula');
    cb = colorbar(ax); cb.Label.String = 'Confidence';
    xlabel(ax, 'X (mm)'); ylabel(ax, 'Y (mm)'); zlabel(ax, 'Z (mm)');
    title(ax, sprintf('Final contacts (n=%d) — colored by confidence', nC));
    axis(ax, 'equal'); grid(ax, 'on'); view(ax, 3);

    % ----- 7) 3D Final by shaft ID + surface -----
    if ~isempty(finalShaftID) && numel(finalShaftID) == nC
        figure('Position', [105 105 600 480], 'Name', [app.PatientIDStr ' — Final by shaft + surface']);
        ax = gca;
        if ~isempty(surfV), patch(ax, 'Vertices', surfV, 'Faces', surfF, 'FaceColor', [0.85 0.85 0.92], 'FaceAlpha', 0.4, 'EdgeColor', 'none'); hold(ax, 'on'); end
        scatter3(ax, finalMM(:,1), finalMM(:,2), finalMM(:,3), sz, finalShaftID, 'filled');
        colormap(ax, 'lines');
        cb = colorbar(ax); cb.Label.String = 'Shaft ID (0=singleton)';
        xlabel(ax, 'X (mm)'); ylabel(ax, 'Y (mm)'); zlabel(ax, 'Z (mm)');
        title(ax, sprintf('Final (n=%d) — colored by shaft', nC));
        axis(ax, 'equal'); grid(ax, 'on'); view(ax, 3);
    end

    % ----- 8) Inter-contact distance distribution -----
    if nC > 1
        D = pdist2(finalMM, finalMM);
        D(logical(eye(size(D)))) = Inf;
        minDist = min(D, [], 2);
        figure('Position', [100 100 500 350], 'Name', [app.PatientIDStr ' — Inter-contact distance']);
        histogram(minDist, 25, 'FaceColor', [0.3 0.6 0.4], 'EdgeColor', 'none');
        xlabel('Min distance to nearest contact (mm)');
        ylabel('Count');
        title(sprintf('Inter-contact spacing (n=%d) — median %.1f mm', nC, median(minDist)));
        grid on;
    end

    % ----- 6) Threshold vs count curve (which threshold “matches” final count) -----
    figure('Position', [120 120 500 350], 'Name', [app.PatientIDStr ' — Threshold curve']);
    plot(ThrHU, NumObj, '-', 'LineWidth', 1.5);
    hold on;
    [~, idxD] = min(abs(NumObj - nC));
    if idxD >= 1 && idxD <= numel(ThrHU)
        plot(ThrHU(idxD), NumObj(idxD), 'ro', 'MarkerSize', 10);
        legend('Raw threshold count', sprintf('Closest to final n=%d', nC), 'Location', 'best');
    end
    xlabel('Threshold (HU)');
    ylabel('# Connected components (raw)');
    title(sprintf('Detection curve vs threshold — final %d contacts in %.1f s', nC, elapsedTime));
    grid on;
end
