function fiberFiltered = ea_filterfiber_roi(ftr, roi, output)
% Filter fibers based on ROI (only keep fibers that traverse the ROI)
%
% Coordinates defined in ftr must be in mm (NOT voxel). It should at least
% have 4 columns in which the 4th column stores the fiber indices.

arguments
    ftr                        % FTR file/structure or fibers coordinates
    roi     {mustBeTextScalar} % ROI path
    output  {mustBeTextScalar} = '' % Optional output file path
end

if isfile(ftr)
    if ~isempty(output)
        ea_mkdir(fileparts(output));
    else
        [~, ftrName] = fileparts(ftr);
        [roiPath, roiName] = ea_niifileparts(GetFullPath(roi));
        output = fullfile(fileparts(roiPath), [ftrName, '_filtered_by_', roiName, '.mat']);
    end
    fprintf('\nLoading fibers...\n')
    ftr = load(ftr);
    fibers = ftr.fibers;
elseif isstruct(ftr)
    fibers = ftr.fibers;
elseif ismatrix(ftr)
    fibers = ftr;
end

% Initialize fiber connectivity status
fibConn = false(length(unique(fibers(:,4))), 1);

disp('Filtering fibers...');

% Find ROI indices
roiNii = ea_load_nii(roi);
roiInd = find(roiNii.img(:));
[xvox, yvox, zvox] = ind2sub(size(roiNii.img), roiInd);
roiMM = ea_vox2mm([xvox, yvox, zvox], roiNii.mat);

% Trim input fibers based on the bounding box
filter = all(fibers(:,1:3)>=min(roiMM),2) & all(fibers(:,1:3)<=max(roiMM), 2);

% Skip further calculation in case ROI is totally not connected
if ~any(filter)
    disp('No fibers found, skipping...');
    fiberFiltered = [];
    return;
end

trimmedFiber = fibers(filter,:);

% Map fibers into ROI voxel space
[trimmedFiberInd, ~, trimmedFiberID] = unique(trimmedFiber(:,4), 'stable');
fibVoxInd = splitapply(@(fib) {ea_mm2uniqueVoxInd(fib, roiNii)}, trimmedFiber(:,1:3), trimmedFiberID);

% Remove outliers
fibVoxInd(cellfun(@(x) any(isnan(x)), fibVoxInd)) = [];
trimmedFiberInd(cellfun(@(x) any(isnan(x)), fibVoxInd)) = [];

% Find connected fibers
connected = cellfun(@(fib) any(ismember(fib, roiInd)), fibVoxInd);

% Binary fiber connection matrix for the T-test method
fibConn(trimmedFiberInd(connected))=1;

% Verbose
ea_cprintf('Comment', '%d out of %d fibers found...\n', sum(fibConn), length(fibConn));

if ~any(connected)
    disp('No fibers found, skipping...');
    fiberFiltered = [];
    return;
end

% Final filtered fibers
fibers = fibers(ismember(fibers(:,4), find(fibConn)), :);

% Adapt fiber indices
oldInd = unique(fibers(:,4), 'stable');
for i=1:length(oldInd)
    fibers(fibers(:,4)==oldInd(i),4) = i;
end

% Construct fiber idx
[~, ~, idx] = unique(fibers(:,4), 'stable');
idx = accumarray(idx,1);

% Handle return value and output
if isstruct(ftr)
    fiberFiltered = ftr;
    fiberFiltered.fibers = fibers;
    fiberFiltered.idx = idx;
    if isfield(fiberFiltered, 'vals')
        fiberFiltered.vals = fiberFiltered.vals(fibConn);
    end

    if ~isempty(output)
        save(output, '-struct', 'fiberFiltered', '-v7.3');
        ea_cprintf('Comment', 'Filtered fibers saved to:\n%s\n', output);
    end
elseif ismatrix(ftr)
    fiberFiltered = fibers;

    if ~isempty(output)
        ftr = struct;    
        ftr.ea_fibformat = '1.1';
        ftr.fibers = fibers;
        ftr.fourindex = 1;
        ftr.idx = idx;
        save(output, '-struct', 'ftr', '-v7.3');
        ea_cprintf('*Comment', 'Filtered fibers saved to:\n%s\n', output);
    end
end

fprintf('Finished!\n');
