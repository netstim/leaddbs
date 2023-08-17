function ea_genvat_wrapper(patientFolder, stimLabel, side, cathode, stimAmplitude, anode, stimType, estimateInTemplate)
% Wrapper to create VTA without GUI.
% Only support single source evenly splitted SimBio model for now.
% Contacts are indexed as 0-N. Case is indexed as -1.

if endsWith(patientFolder, filesep)
    patientFolder = patientFolder(1:end-1);
end

[~, patientName] = fileparts(patientFolder);

% Create native space output folder
nativeStimFolder = fullfile(patientFolder, 'stimulations', ea_nt(1), stimLabel);
if ~isfolder(nativeStimFolder)
    mkdir(nativeStimFolder);
end

% Create template space output folder
templateStimFolder = fullfile(patientFolder, 'stimulations', ea_nt(0), stimLabel);
if ~isfolder(templateStimFolder)
    mkdir(templateStimFolder);
end

% Assume voltage stimulation by default
if ~exist('anode', 'var') || isempty(anode) || strcmpi(anode, 'case')
    anode = -1; % Case is anode
end

% Assume voltage stimulation by default
if ~exist('stimType', 'var') || isempty(stimType)
    stimType = 'voltage';
end

% Estimate in native space by default
if ~exist('estimateInTemplate', 'var') || isempty(estimateInTemplate)
    estimateInTemplate = 0;
end

% Load template stimparameter structure
load([ea_getearoot, 'common', filesep, 'stimparameters_template.mat'], 'S');

% Get patient options
options = ea_getptopts(patientFolder);

% Save original native flag
options.orignative = options.native;

% Set stimulation label
S.label = stimLabel;

% Set native flag depending on estimation space
if ~estimateInTemplate
    options.native = 1;
end

switch lower(side)
    case {0, 'r', 'right'}
        sideChar = 'R';
        sideInd = 1;
        sideOffSet = 0; % k0-k7, no index offset
    case {1, 'l', 'left'}
        sideChar = 'L';
        sideInd = 2;
        sideOffSet = 8; % k8-k15, offset the indices from 0-7 to 8-15
end

switch lower(stimType)
    case {1, 'voltage'}
        stimType = 1;
    case {2, 'current'}
        stimType = 2;
end

% Set stimulation amplitude and type
S.amplitude{sideInd}(1) = stimAmplitude;
eval(['S.',sideChar,'s1.amp = stimAmplitude;']);
eval(['S.',sideChar,'s1.va = stimType;']);

% Set active contact
S.activecontacts{sideInd}(setdiff(union(anode, cathode),-1)+1) = 1;

% Set anode parameters
for i=1:numel(anode)
    if anode(i) == -1
        eval(['S.',sideChar,'s1.case.perc = 1/numel(anode)*100;']);
        eval(['S.',sideChar,'s1.case.pol = 2;']);
    else
        eval(['S.',sideChar,'s1.k',num2str(anode(i)+sideOffSet),'.perc = 1/numel(anode)*100;']);
        eval(['S.',sideChar,'s1.k',num2str(anode(i)+sideOffSet),'.pol = 2;']);
    end
end

% Set cathode parameters
for i=1:numel(cathode)
    if cathode(i) == -1
        eval(['S.',sideChar,'s1.case.perc = 1/numel(cathode)*100;']);
        eval(['S.',sideChar,'s1.case.pol = 1;']);
    else
        eval(['S.',sideChar,'s1.k',num2str(cathode(i)+sideOffSet),'.perc = 1/numel(cathode)*100;']);
        eval(['S.',sideChar,'s1.k',num2str(cathode(i)+sideOffSet),'.pol = 1;']);
    end
end

% Save stimulation parameters in template space output folder
S.template = 'direct';
save(fullfile(templateStimFolder, [patientName, '_desc-stimparameters.mat']), 'S');

% Save stimulation parameters in native space output folder
S.template = 'warp';
save(fullfile(nativeStimFolder, [patientName, '_desc-stimparameters.mat']), 'S');

% Fix missing atlasset in options
options.atlasset = 'DISTAL Nano (Ewert 2017)';

% Run VTA calculation
ea_genvat_horn([], S, sideInd, options, stimLabel);
