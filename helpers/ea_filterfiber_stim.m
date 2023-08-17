function fiberFiltered = ea_filterfiber_stim(ftr, coords, S, type, factor, ref)
% Filter fibers based on the active contacts and stimulation amplitudes

fprintf('\nCollecting stimulation parameters...\n')

% Active contacts indices
if iscell(S) % stimSetMode, stimProtocol (cell of csv files) provided
    stimProtocol = S;
    stimProtocol = cellfun(@(f) table2array(readtable(f,'NumHeaderLines',1)), stimProtocol, 'Uni', 0)';
    % special case for unilateral StimSets
    % assign a null stim protocol for the other side
    if size(stimProtocol,2) == 1
        disp("StimProtocol exists only for one side")
        N_contacts = size(stimProtocol{1,1},2);
        [~,SetName,~] = fileparts(S);
        if strcmp(SetName, 'Current_protocols_0')
            N_contacts = size(stimProtocol{1,1},2);
            stimProtocol{1,2} = zeros(1,N_contacts);
        elseif strcmp(SetName, 'Current_protocols_1')
            % swap sides
            stimProtocol{1,2} = stimProtocol{1,1};
            stimProtocol{1,1} = zeros(1,N_contacts);
        else
            disp("Unrecongnized protocol file name, exitting")
            return
        end
    end


    activeContacts = cell(size(stimProtocol));
    for i=1:length(activeContacts)
        activeContacts{i} = find(~isnan(max(stimProtocol{i}))); % Find contact with stimulation input
    end

    
    stimAmplitudes = cell(size(stimProtocol));
    for i=1:length(stimAmplitudes)
        stimAmplitudes{i} = max(abs(stimProtocol{i})); % Use maximum absolute amplitude
    end
    
else % normal mode
    if ischar(S) && isfile(S) % stimparameters.mat provided
        load(S, 'S');
    end

    activeContacts = cell(size(S.activecontacts));
    for i=1:length(activeContacts)
        activeContacts{i} = find(S.activecontacts{i});
    end
    
    % define the stim vector as in OSS-DBS (sources are merged)
    stimAmplitudes = cell(size(S.amplitude));
    eleNum = length(coords); % Number of electrodes
    conNum = cellfun(@(x) size(x,1), coords); % Number of contacts per electrode
    conNum = conNum(find(conNum, 1));
    stimVector = ea_getStimVector(S, eleNum, conNum);
    for side = 1:size(S.amplitude,2)
        for cnt = 1:size(stimVector(side,:),2)
            if isnan(stimVector(side,cnt))
                stimAmplitudes{side}(cnt) = 0.0;
            else 
                stimAmplitudes{side}(cnt) = abs(stimVector(side,cnt)); % sign does not matter for Kuncel-VTA
            end
        end
    end
end

% Active contacts coordinates
stimCoords = cell(size(coords));
for i=1:length(stimCoords)
    stimCoords{i} = coords{i}(activeContacts{i},:);
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

load([ea_space, 'spacedef.mat'], 'spacedef');
primarytemplate = spacedef.templates{1};

% Reference image when construct the spherical ROI, use MNI t1 by default
if ~exist('ref','var')
    ref = [ea_space, primarytemplate, '.nii'];
end

% Check if fibers pass through the ROI
fibConn = zeros(length(ftr.idx), length(stimAmplitudes));
fiberFiltered = cell(size(stimAmplitudes));
for i=1:length(radius)
    if ~radius{i}
        disp('No stimulation found, skipping...');
    else
        fprintf('\nConstructing spherical ROI...\n');
        sphereROI = ea_spherical_roi([],stimCoords{i}, radius{i}, 0, ref);

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
    U(U<0.1) = 0.1;  % fix for algorithm-based protocols
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
