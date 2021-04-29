function ea_split_discfibers_by_vta(vatlist, tractMATFile, outputFolder)
% Split fibers in the discfiber atlas into TRACT_[In|Out]_VTA files
% depending on whether the fibers are connected to the VTAs or not.
%
% TRACT_In_VTA contains fibers which are connected to the VTA.
% TRACT_Out_VTA contains fibers which are not connected to the VTA.
%
% Resulted files can be visualized in the 3D visualization window from the
% menu: "Add Objects" -> "Open Tracts".

% List of VTA nifti files
if ischar(vatlist)
    vatlist = {vatlist};
end

% Output to current path by default
if ~exist('outputFolder', 'var')
    outputFolder = pwd;
else
    ea_mkdir(outputFolder);
end

% Load tracts from discfiber atlas mat file
disp('Load discfiber atlas...');
load(tractMATFile, 'fibcell');
[~, tractName] = fileparts(tractMATFile);

% Stack the fibers and vals from the chosen side(s) all together
fibcell = vertcat(fibcell{:});

% Reform fibers to N*4 matrix
fibers = vertcat(fibcell{:});
idx = cellfun(@(fib) size(fib,1), fibcell);
fibers(:,4) = repelem(1:numel(fibcell), idx);

% Threshold for E-field
prefs = ea_prefs;
thresh = prefs.machine.vatsettings.horn_ethresh*1000;

numVAT = numel(vatlist);

% Define fiber connection matrix
fibConnBin = zeros(numel(fibcell),1);

ea_fibformat = '1.0';
fourindex = 1;

for v = 1:numVAT
    disp(['VAT ', num2str(v, ['%0',num2str(numel(num2str(numVAT))),'d']), '/', num2str(numVAT), '...']);
    vat = ea_load_nii(vatlist{v});
    [~, vtaName] = ea_niifileparts(vatlist{v});

    % Threshold the vat efield
    if numel(unique(vat.img(:))) ~= 2
        vatInd = find(vat.img(:)>thresh);
    else
        vatInd = find(vat.img(:));
    end

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
    fibVoxInd = splitapply(@(fib) {ea_mm2uniqueVoxInd(fib, vat)}, trimmedFiber(:,1:3), trimmedFiberID);

    % Remove outliers
    fibVoxInd(cellfun(@(x) any(isnan(x)), fibVoxInd)) = [];
    trimmedFiberInd(cellfun(@(x) any(isnan(x)), fibVoxInd)) = [];

    % Find connected fibers
    connected = cellfun(@(fib) any(ismember(fib, vatInd)), fibVoxInd);

    % Binary fiber connection matrix
    fibConnBin(trimmedFiberInd(connected))=1;
    fibConnBin = logical(fibConnBin);

    % Get connected fibers
    if sum(fibConnBin)
        disp([num2str(sum(fibConnBin)), ' fibers connected to ', vtaName, '!']);
        fibcellIn = fibcell(fibConnBin);
        % Reform fibers to N*4 matrix
        fibers = vertcat(fibcellIn{:});
        idx = cellfun(@(fib) size(fib,1), fibcellIn);
        fibers(:,4) = repelem(1:numel(fibcellIn), idx);

        % Save fibers
        save([outputFolder, filesep, tractName, '_In_', vtaName], 'ea_fibformat', 'fourindex', 'fibers', 'idx');
    else
        disp(['0 fibers connected to ', vtaName, '!']);
    end

    % Get unconnected fibers
    if sum(~fibConnBin)
        disp([num2str(sum(~fibConnBin)), ' fibers not connected to ', vtaName, '!']);
        fibcellOut = fibcell(~fibConnBin);
        % Reform fibers to N*4 matrix
        fibers = vertcat(fibcellOut{:});
        idx = cellfun(@(fib) size(fib,1), fibcellOut);
        fibers(:,4) = repelem(1:numel(fibcellOut), idx);

        % Save fibers
        save([outputFolder, filesep, tractName, '_Out_', vtaName], 'ea_fibformat', 'fourindex', 'fibers', 'idx');
    else
        disp(['0 fibers unconnected to ', vtaName, '!']);
    end
end
