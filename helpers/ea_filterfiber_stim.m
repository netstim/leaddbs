function fiberFiltered = ea_filterfiber_stim(ftr, coords, S, type, factor)
% Filter fibers based on the active contacts and stimulation amplitudes

% Load stimparameters in case needed
if ischar(S) && isfile(S)
    load(S, 'S');
end

fprintf('\nCollecting stimulation parameters...\n')

% Active contacts indices
activeContacts = cell(size(S.activecontacts));
for i=1:length(activeContacts)
    activeContacts{i} = find(S.activecontacts{i});
end

% Active contacts coordinates
stimCoords = cell(size(coords));
for i=1:length(stimCoords)
    stimCoords{i} = coords{i}(activeContacts{i},:);
end

% Stimulation amplitude, V or mA.
stimAmplitudes = cell(size(S.amplitude));
for i=1:length(stimAmplitudes)
    stimAmplitudes{i} = zeros(1, size(coords{1},1));
end

for side = 1:length(S.amplitude)
    % Set stimulation source label and contacts labels
    switch side
        case 1
            sideCode = 'R';
            cntlabel = {'k0','k1','k2','k3','k4','k5','k6','k7'};
        case 2
            sideCode = 'L';
            cntlabel = {'k8','k9','k10','k11','k12','k13','k14','k15'};
    end

    % Index of stimulation source, only support 1 input source for now
    sourceIndex = find(S.amplitude{side},1);
    if ~isempty(sourceIndex)
        stimSource = S.([sideCode, 's', num2str(sourceIndex)]);

        % Collect stimulation amplitudes
        for cnt = 1:length(stimAmplitudes{side})
            if S.activecontacts{side}(cnt)
                stimAmplitudes{side}(cnt) = S.amplitude{side}(sourceIndex)*stimSource.(cntlabel{cnt}).perc/100;
            end
        end
    end
end

% Only keep amplitudes for active contacts
for i=1:length(stimAmplitudes)
    stimAmplitudes{i} = stimAmplitudes{i}(:, activeContacts{i});
end

% Define anonymous function to calc radius
if ~exist('type', 'var')
    type = 'kuncel';
end

switch lower(type)
    case 'kuncel'
        calcr = @(U) kuncel08_eq1(U);
        fprintf('\nEstimating radius based on Kuncel et al. 2008...\n');
    case 'maedler'
        calcr = @(U) maedler12_eq3(U);
        fprintf('\nEstimating radius based on Maedler et al. 2012...\n')
end

% Calculate radius
if ~exist('factor', 'var')
    factor = 1;  % Multiply the radius from the estimation by a fixed factor
end
radius = cell(size(stimAmplitudes));
for i=1:length(radius)
    radius{i} = calcr(stimAmplitudes{i});
    disp(['Radius (mm): ', num2str(radius{i}), ' x ', num2str(factor)]);
    % Triple the radius since the model tends to underestimate activation
    radius{i} = radius{i} * factor;
end

% Load fiber connectome
if ischar(ftr) && isfile(ftr)
    fprintf('\nLoading fibers...\n');
    ftr = load(ftr);
end

% Convert voxel coordinates to mm in case needed
if isfield(ftr, 'voxmm') && strcmp(ftr.voxmm, 'vox')
    fprintf('\nConverting voxel coordinates to mm...\n')
    if isfield(ftr, 'mat')
        ftr.fibers(:,1:3) = ea_vox2mm(ftr.fibers(:,1:3), ftr.mat);
    else
        error('Fibers stored in voxel space but affine matrix doesn''t exist!');
    end
end

% Check if fibers pass through the ROI
fibConn = zeros(length(ftr.idx), length(stimAmplitudes));
fiberFiltered = cell(size(stimAmplitudes));
for i=1:length(radius)
    if ~radius{i}
        disp('No stimulation found, skipping...');
    else
        fprintf('\nConstructing spherical ROI...\n');
        sphereROI = ea_spherical_roi([],stimCoords{i}, radius{i}, 0);

        % Find indices within sphere ROI region
        sphereROIInd = find(sphereROI.img(:));

        disp('Filtering fibers...');

        % Trim connectome fibers
        [xvox, yvox, zvox] = ind2sub(size(sphereROI.img), sphereROIInd);
        sphereROImm = ea_vox2mm([xvox, yvox, zvox], sphereROI.mat);
        filter = all(ftr.fibers(:,1:3)>=min(sphereROImm),2) & all(ftr.fibers(:,1:3)<=max(sphereROImm), 2);

        % Skip further calculation in case VAT is totally not connected
        if ~any(filter)
            disp('No fibers found, skipping...');
            continue;
        end

        trimmedFiber = ftr.fibers(filter,:);

        % Map mm connectome fibers into VAT voxel space
        [trimmedFiberInd, ~, trimmedFiberID] = unique(trimmedFiber(:,4), 'stable');
        fibVoxInd = splitapply(@(fib) {ea_mm2uniqueVoxInd(fib, sphereROI)}, trimmedFiber(:,1:3), trimmedFiberID);

        % Remove outliers
        fibVoxInd(cellfun(@(x) any(isnan(x)), fibVoxInd)) = [];
        trimmedFiberInd(cellfun(@(x) any(isnan(x)), fibVoxInd)) = [];

        % Find connected fibers
        connected = cellfun(@(fib) any(ismember(fib, sphereROIInd)), fibVoxInd);

        % Binary fiber connection matrix for the T-test method
        fibConn(trimmedFiberInd(connected), i)=1;

        fprintf('%d out of %d fibers found...\n', sum(fibConn(:,i)), length(ftr.idx));
        fiberFiltered{i}.fibers = ftr.fibers(ismember(ftr.fibers(:,4), find(fibConn(:,i))), :);
        fiberFiltered{i}.idx = ftr.idx(logical(fibConn(:,i)));
    end
end

fprintf('\nFinished!\n\n');


function r = kuncel08_eq1(U)
% This function calculates the  radius of Volume of Activated Tissue for
% stimulation settings U (Kuncel 2008).

r = 0;
if U
    k = 0.22;
    Uo = 0.1;
    r = sqrt((U-Uo)/k);
end


function r = maedler12_eq3(U)
% This function calculates the  radius of Volume of Activated Tissue for
% stimulation settings U (Maedler 2012). Clinical measurements of DBS
% electrode impedance typically range from 500-1500 Ohm (Butson 2006).

Im = 1000;

r = 0;
if U
    k1 =-1.0473;
    k3 = 0.2786;
    k4 = 0.0009856;
    r = -(k4*Im - sqrt(k4^2*Im^2 + 2*k1*k4*Im + k1^2 + 4*k3*U) + k1)/(2*k3);
end
