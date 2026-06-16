function out = LeG_autoNameElecs(app)
% LeG_autoNameElecs  Group detected contacts into shafts and assign labels.

coords = pickCoordsAutoName(app, {'ElecXYZRaw', 'ElecXYZProjRaw'});
dispCoords = pickCoordsAutoName(app, {'ElecXYZProjRaw', 'ElecXYZRaw'});
nContacts = size(coords, 1);
if nContacts == 0
    error('LeG_autoNameElecs:NoContacts', 'No electrodes are available to name.');
end
if size(dispCoords, 1) ~= nContacts
    dispCoords = coords;
end

spacingMm = estimateSpacingAutoName(coords);
shafts = detectShaftsAutoName(coords, spacingMm);
shafts = attachRemainingContactsAutoName(coords, shafts, spacingMm);
shafts = addSingletonShaftsAutoName(coords, shafts);
%shafts = finalizeShaftsAutoName(coords, dispCoords, shafts, getSurfaceVerticesAutoName(app));
%EG added this
surfaceVerts = getSurfaceVerticesAutoName(app);

shafts = finalizeShaftsAutoName(coords, dispCoords, shafts, surfaceVerts);

shafts = forceAllShaftsDistalToProximalAutoName(coords, shafts, surfaceVerts);
%end
nShafts = numel(shafts);
counts = arrayfun(@(x) numel(x.indices), shafts);
maxContacts = max(counts);

elecMapRaw = cell(nContacts, 3);
elecMapRaw(:, 1) = {'NaN'};
elecMapRaw(:, 2:3) = num2cell(nan(nContacts, 2));
labelMap = repmat({'NaN'}, maxContacts, nShafts);
channelMap1 = nan(maxContacts, nShafts);
channelMap2 = nan(maxContacts, nShafts);
shaftMembership = nan(nContacts, 1);

chan = 1;
leadNames = cell(nShafts, 1);
for s = 1:nShafts
    leadNames{s} = shaftNameAutoName(s);
    idx = shafts(s).indices(:);
    rows = maxContacts:-1:(maxContacts-numel(idx)+1);
    for c = 1:numel(idx)
        label = sprintf('%s%d', leadNames{s}, c);
        elecMapRaw(idx(c), :) = {label, chan, nan};
        labelMap{rows(c), s} = label;
        channelMap1(rows(c), s) = chan;
        shaftMembership(idx(c)) = s;
        chan = chan + 1;
    end
end

out = struct();
out.ElecMapRaw = elecMapRaw;
out.LabelMap = labelMap;
out.ChannelMap1 = channelMap1;
out.ChannelMap2 = channelMap2;
out.DepthElecRaw = (1:nContacts)';
out.MicroElecRaw = [];
out.RefElecRaw = [];
out.GndElecRaw = [];
out.ShaftMembership = shaftMembership;
out.ShaftNames = leadNames;
out.SpacingMm = spacingMm;
out.Shafts = shafts;
end

function coords = pickCoordsAutoName(app, fieldOrder)
coords = [];
for k = 1:numel(fieldOrder)
    fld = fieldOrder{k};
    if (isobject(app) && isprop(app, fld)) || (isstruct(app) && isfield(app, fld))
        val = app.(fld);
        if isnumeric(val) && size(val, 2) == 3 && ~isempty(val)
            coords = double(val);
            break;
        end
    end
end
if isempty(coords)
    coords = zeros(0, 3);
end
end

function verts = getSurfaceVerticesAutoName(app)

verts = [];

surfaceFields = {'ProjSurfRaw', 'BrainSurfRaw'};

for k = 1:numel(surfaceFields)
    fld = surfaceFields{k};

    if (isobject(app) && isprop(app, fld)) || (isstruct(app) && isfield(app, fld))
        surf = app.(fld);

        if isstruct(surf) && isfield(surf, 'vertices') && size(surf.vertices, 2) == 3
            verts = double(surf.vertices);
            return;
        end
    end
end

end

function spacingMm = estimateSpacingAutoName(coords)
n = size(coords, 1);
if n < 2
    spacingMm = 5.0;
    return;
end

nn = inf(n, 1);
for i = 1:n
    d = sqrt(sum((coords - coords(i, :)).^2, 2));
    d(i) = inf;
    nn(i) = min(d);
end

nn = nn(isfinite(nn) & nn > 0);
if isempty(nn)
    spacingMm = 5.0;
else
    spacingMm = median(nn);
end
spacingMm = min(max(spacingMm, 2.0), 10.0);
end

function shafts = detectShaftsAutoName(coords, spacingMm)
shafts = emptyShaftsAutoName();
remaining = (1:size(coords, 1))';

while numel(remaining) >= 2
    localPts = coords(remaining, :);
    pairs = candidatePairsAutoName(localPts, spacingMm);
    if isempty(pairs)
        break;
    end

    bestIdx = [];
    bestScore = -inf;
    for p = 1:size(pairs, 1)
        [idx, score] = evaluatePairAutoName(localPts, pairs(p, 1), pairs(p, 2), spacingMm);
        if score > bestScore
            bestScore = score;
            bestIdx = idx;
        end
    end

    if numel(bestIdx) < 2 || ~isfinite(bestScore)
        break;
    end

    globalIdx = remaining(bestIdx);
    shafts(end+1) = refitShaftAutoName(coords, globalIdx); %#ok<AGROW>

    keep = true(numel(remaining), 1);
    keep(bestIdx) = false;
    remaining = remaining(keep);
end
end

function pairs = candidatePairsAutoName(points, spacingMm)
n = size(points, 1);
pairs = zeros(0, 2);
if n < 2
    return;
end

pairDistMax = max(18.0, 4.0 * spacingMm);
kPairs = min(4, n-1);
D = inf(n, n);
for i = 1:n-1
    d = sqrt(sum((points(i+1:end, :) - points(i, :)).^2, 2));
    D(i, i+1:end) = d;
    D(i+1:end, i) = d;
end

mask = false(n, n);
for i = 1:n
    [sd, ord] = sort(D(i, :), 'ascend');
    ord = ord(sd <= pairDistMax);
    if isempty(ord)
        continue;
    end
    ord = ord(1:min(kPairs, numel(ord)));
    mask(i, ord) = true;
end

[ii, jj] = find(triu(mask | mask', 1));
pairs = [ii, jj];
if isempty(pairs)
    [ii, jj] = find(triu(true(n), 1));
    pairs = [ii, jj];
end
end

function [bestIdx, bestScore] = evaluatePairAutoName(points, i1, i2, spacingMm)
bestIdx = [];
bestScore = -inf;

dir = points(i2, :) - points(i1, :);
ndir = norm(dir);
if ndir < eps
    return;
end
dir = dir / ndir;

t = (points - points(i1, :)) * dir';
proj = points(i1, :) + t * dir;
perp = sqrt(sum((points - proj).^2, 2));

inlierRadius = min(2.5, max(1.25, 0.35 * spacingMm + 0.2));
gapThr = max(8.5, 2.25 * spacingMm);
inliers = find(perp <= inlierRadius);
if numel(inliers) < 2
    return;
end

[tSorted, ord] = sort(t(inliers), 'ascend');
inliers = inliers(ord);
cutIdx = [0; find(diff(tSorted) > gapThr); numel(inliers)];
for k = 1:numel(cutIdx)-1
    idx = inliers(cutIdx(k)+1:cutIdx(k+1));
    if numel(idx) < 2
        continue;
    end
    span = t(idx(end)) - t(idx(1));
    if numel(idx) == 2 && span < max(4.0, 0.85 * spacingMm)
        continue;
    end

    score = 10 * numel(idx) + span / max(spacingMm, 1e-3) - ...
        2.0 * mean(perp(idx)) / max(inlierRadius, 1e-3);
    if score > bestScore
        bestScore = score;
        bestIdx = idx(:);
    end
end
end

function shafts = attachRemainingContactsAutoName(coords, shafts, spacingMm)
if isempty(shafts)
    return;
end

attachRadius = min(2.8, max(1.4, 0.45 * spacingMm));
spanMargin = max(5.0, 1.4 * spacingMm);
used = false(size(coords, 1), 1);
for k = 1:numel(shafts)
    used(shafts(k).indices) = true;
end

remaining = find(~used);
for i = 1:numel(remaining)
    idx = remaining(i);
    pt = coords(idx, :);
    bestShaft = 0;
    bestPerp = inf;
    for s = 1:numel(shafts)
        [perpDist, tVal] = shaftDistanceAutoName(pt, shafts(s));
        if perpDist > attachRadius
            continue;
        end
        if tVal < shafts(s).tRange(1) - spanMargin || tVal > shafts(s).tRange(2) + spanMargin
            continue;
        end
        if perpDist < bestPerp
            bestPerp = perpDist;
            bestShaft = s;
        end
    end
    if bestShaft > 0
        shafts(bestShaft) = refitShaftAutoName(coords, [shafts(bestShaft).indices(:); idx]);
    end
end
end

function shafts = addSingletonShaftsAutoName(coords, shafts)
used = false(size(coords, 1), 1);
for k = 1:numel(shafts)
    used(shafts(k).indices) = true;
end

remaining = find(~used);
for i = 1:numel(remaining)
    shafts(end+1) = refitShaftAutoName(coords, remaining(i)); %#ok<AGROW>
end
end

function shafts = finalizeShaftsAutoName(coords, dispCoords, shafts, surfaceVerts)
for k = 1:numel(shafts)
    shafts(k) = refitShaftAutoName(coords, shafts(k).indices);
    localIdx = orientShaftAutoName(coords(shafts(k).indices, :), surfaceVerts);
    shafts(k).indices = shafts(k).indices(localIdx);
    shaftPts = dispCoords(shafts(k).indices, :);
    ctr = mean(shaftPts, 1);
    shafts(k).sortKey = [ctr(1), -ctr(2), -ctr(3), -numel(shafts(k).indices)];
end

if numel(shafts) > 1
    [~, ord] = sortrows(cat(1, shafts.sortKey));
    shafts = shafts(ord);
end
end

function order = orientShaftAutoName(points, surfaceVerts)
[~, ~, tSorted, order] = shaftAxisAutoName(points);
points = points(order, :);
if size(points, 1) < 2
    return;
end

flipFlag = false;
if ~isempty(surfaceVerts)
    d1 = pointSurfaceDistanceAutoName(points(1, :), surfaceVerts);
    d2 = pointSurfaceDistanceAutoName(points(end, :), surfaceVerts);
    if d1 < d2 - 0.5
        flipFlag = true;
    elseif abs(d1 - d2) <= 0.5
        flipFlag = norm(points(1, :)) > norm(points(end, :));
    end
else
    flipFlag = norm(points(1, :)) > norm(points(end, :));
end

if flipFlag
    order = flipud(order);
end
end

function d = pointSurfaceDistanceAutoName(point, surfaceVerts)
d = sqrt(min(sum((surfaceVerts - point).^2, 2)));
end

function [perpDist, tVal] = shaftDistanceAutoName(point, shaft)
if numel(shaft.indices) == 1
    perpDist = norm(point - shaft.mu);
    tVal = 0;
    return;
end

dp = point - shaft.mu;
tVal = dp * shaft.dir';
proj = shaft.mu + tVal * shaft.dir;
perpDist = norm(point - proj);
end

function shaft = refitShaftAutoName(coords, idx)
idx = unique(idx(:), 'stable');
pts = coords(idx, :);
[dir, mu, tSorted, ord, span] = shaftAxisAutoName(pts);

shaft = struct();
shaft.indices = idx(ord);
shaft.dir = dir;
shaft.mu = mu;
shaft.tRange = [tSorted(1), tSorted(end)];
shaft.span = span;
shaft.sortKey = zeros(1, 4);
end

function [dir, mu, tSorted, ord, span] = shaftAxisAutoName(points)
mu = mean(points, 1);
if size(points, 1) == 1
    dir = [1, 0, 0];
    tSorted = 0;
    ord = 1;
    span = 0;
    return;
end

X = points - mu;
[~, ~, V] = svd(X, 0);
dir = V(:, 1)';
t = X * dir';
[tSorted, ord] = sort(t, 'ascend');
span = tSorted(end) - tSorted(1);
end

function name = shaftNameAutoName(idx)
letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ';
name = '';
while idx > 0
    remIdx = mod(idx - 1, 26) + 1;
    name = [letters(remIdx), name]; %#ok<AGROW>
    idx = floor((idx - 1) / 26);
end
end

function shafts = emptyShaftsAutoName()
shafts = struct('indices', {}, 'dir', {}, 'mu', {}, 'tRange', {}, ...
    'span', {}, 'sortKey', {});
end
%EG added this function
function shafts = forceAllShaftsDistalToProximalAutoName(coords, shafts, surfaceVerts)

for s = 1:numel(shafts)

    idx = shafts(s).indices(:);

    if numel(idx) < 2
        continue;
    end

    p1 = coords(idx(1), :);
    p2 = coords(idx(end), :);

    d1 = norm(p1);  % distance of first contact from [0 0 0]
    d2 = norm(p2);  % distance of last contact from [0 0 0]

    % Distal should be closer to MNI origin.
    % If the first contact is farther from origin than the last contact,
    % flip this shaft only.
    if d1 > d2
        shafts(s).indices = flipud(idx);
    end

end

end