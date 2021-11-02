function [vatFibScoreBin, vatFibScoreSum, vatFibScoreMean, vatFibScorePeak, vatFibScore5Peak] = ea_discfibers_vtascore(vatlist, atlas, side, posneg, threshold)
% Calculate VAT fiber connection scores based on the discfiber atlas

if ischar(vatlist)
    vatlist = {vatlist};
end

disp('Load discfiber atlas...');
load(atlas, 'fibcell', 'vals');

if ~exist('side', 'var')
    side = 'both';
end

if ~exist('posneg', 'var')
    posneg = 'both';
end

% Predict using fiberss from only one side
switch lower(side)
    case 'right'
        if isempty(fibcell{1})
            error('No fiber exists for the right side!')
        else
            fibcell = fibcell(1);
            vals = vals(1);
        end
    case 'left'
        if numel(fibcell)==1 || isempty(fibcell{2})
            error('No fiber exists for the right side!')
        else
            fibcell = fibcell(2);
            vals = vals(2);
        end
end

% Stack the fibers and vals from the chosen side(s) all together
fibcell = vertcat(fibcell{:});
vals = vertcat(vals{:});

% Predict using only positive/negative fibers
switch lower(posneg)
    case {'pos', 'positive'}
        if all(vals<0)
            error('No positive fiber exists!')
        else
            fibcell = fibcell(vals>0);
            vals = vals(vals>0);
        end
    case {'neg', 'negative'}
        if all(vals>0)
            error('No negative fiber exists!')
        else
            fibcell = fibcell(vals<0);
            vals = vals(vals<0);
        end
end

% Reform fibers to N*4 matrix
fibers = vertcat(fibcell{:});
idx = cellfun(@(fib) size(fib,1), fibcell);
fibers(:,4) = repelem(1:numel(fibcell), idx);

% Use E-field threshold by default. Otherwise, the input vatlist should be
% properly binarized.
if ~exist('threshold', 'var')
    prefs = ea_prefs;
    threshold = prefs.machine.vatsettings.horn_ethresh*1000;
end

numVAT = numel(vatlist);

% Define fiber connection matrix
fibConnBin = zeros(numel(fibcell), numVAT);
fibConnSum = zeros(numel(fibcell), numVAT);
fibConnMean = zeros(numel(fibcell), numVAT);
fibConnPeak = zeros(numel(fibcell), numVAT);
fibConn5Peak = zeros(numel(fibcell), numVAT);

for pt = 1:numVAT
    disp(['VAT ', num2str(pt, ['%0',num2str(numel(num2str(numVAT))),'d']), '/', num2str(numVAT), '...']);
    vat = ea_load_nii(vatlist{pt});

    % Threshold the vat efield
    if numel(unique(vat.img(:))) ~= 2
        vatInd = find(vat.img(:)>threshold);
    else
        vatInd = find(vat.img(:));
    end

    if isempty(vatInd)
        warning('off', 'backtrace');
        warning('Skip empty VTA %s ...', vat.fname);
        warning('on', 'backtrace');
        continue;
    end

    % Trim connectome fibers
    [xvox, yvox, zvox] = ind2sub(size(vat.img), vatInd);
    vatmm = ea_vox2mm([xvox, yvox, zvox], vat.mat);
    filter = all(fibers(:,1:3)>=min(vatmm),2) & all(fibers(:,1:3)<=max(vatmm), 2);

    % Skip further calculation in case VAT is totally not connected
    if ~any(filter)
        warning('off', 'backtrace');
        warning('Skip unconnected VTA %s ...', vat.fname);
        warning('on', 'backtrace');
        continue;
    end

    trimmedFiber = fibers(filter,:);

    % Map mm connectome fibers into VAT voxel space
    [trimmedFiberInd, ~, trimmedFiberID] = unique(trimmedFiber(:,4), 'stable');
    fibVoxInd = splitapply(@(fib) {ea_mm2uniqueVoxInd(fib, vat)}, trimmedFiber(:,1:3), trimmedFiberID);

    % Remove outliers
    fibVoxInd(cellfun(@(x) any(isnan(x)), fibVoxInd)) = [];
    trimmedFiberInd(cellfun(@(x) any(isnan(x)), fibVoxInd)) = [];

    % Find connected fibers
    connected = cellfun(@(fib) any(ismember(fib, vatInd)), fibVoxInd);

    % Binary fiber connection matrix for the T-test method
    fibConnBin(trimmedFiberInd(connected), pt)=1;

    % Checck intersection between vat and the connected fibers
    intersection = cellfun(@(fib) vat.img(intersect(fib, vatInd)), fibVoxInd(connected), 'Uni', 0);

    % Fiber connection matrix for the Spearman's correlation method
    fibConnSum(trimmedFiberInd(connected), pt) = cellfun(@sum, intersection);
    fibConnMean(trimmedFiberInd(connected), pt) = cellfun(@mean, intersection);
    fibConnPeak(trimmedFiberInd(connected), pt) = cellfun(@max, intersection);
    fibConn5Peak(trimmedFiberInd(connected), pt) = cellfun(@(x) mean(maxk(x,ceil(0.05*numel(x)))), intersection);
end

% Calc fiber scores for VATs
vatFibScoreBin = ea_nansum(vals.*fibConnBin)';
vatFibScoreSum = ea_nansum(vals.*fibConnSum)';
vatFibScoreMean = ea_nansum(vals.*fibConnMean)';
vatFibScorePeak = ea_nansum(vals.*fibConnPeak)';
vatFibScore5Peak = ea_nansum(vals.*fibConn5Peak)';
