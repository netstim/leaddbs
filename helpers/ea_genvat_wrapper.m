function ea_genvat_wrapper(patientFolder, stimLabel, side, cathode, stimAmplitude, anode, stimType, estimateInTemplate)
% Wrapper to create VTA without GUI.
% Only support single source evenly splitted SimBio model for now.
% Contacts are indexed as 1-N. Case is indexed as 0.

patientFolder = erase(GetFullPath(patientFolder), filesep + textBoundary("end"));
if ~isfolder(patientFolder)
    ea_error('Patient folder not found!', showdlg=false, simpleStack=true);
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

% Get patient options
options = ea_getptopts(patientFolder);
options.patientname = ['sub-', options.subj.subjId];

% Initialize stimparameter structure
S = ea_initializeS(stimLabel, options);

% Save original native flag
options.orignative = options.native;

% Set native flag depending on estimation space
if ~estimateInTemplate
    options.native = 1;
end

switch lower(side)
    case {0, 'r', 'right'}
        sideChar = 'R';
        sideInd = 1;
    case {1, 'l', 'left'}
        sideChar = 'L';
        sideInd = 2;
end

switch lower(stimType)
    case {1, 'voltage'}
        stimType = 1;
    case {2, 'current'}
        stimType = 2;
end

% Set stimulation amplitude and type
S.amplitude{sideInd}(1) = stimAmplitude;
S.([sideChar,'s1']).amp = stimAmplitude;
S.([sideChar,'s1']).va = stimType;

% Set active contact
S.activecontacts{sideInd}(setdiff(union(anode, cathode),0)) = 1;

% Check anode
if ismember(0, anode) && numel(anode)>1 % When case is anode, it has to be the only anode
    ea_error('When case is anode, it has to be the only anode!', showdlg=false, simpleStack=true);
end

% Check cathode
if ismember(0, cathode) % Case cannot be cathode
    ea_error('Case cannot be cathode!', showdlg=false, simpleStack=true);
end

% Set case
if ismember(0, anode)
    S.([sideChar,'s1']).case.perc = 100;
    S.([sideChar,'s1']).case.pol = 2;
else
    for i=1:numel(anode)
        S.([sideChar,'s1']).(['k',num2str(anode(i))]).perc = 1/numel(anode)*100;
        S.([sideChar,'s1']).(['k',num2str(anode(i))]).pol = 2;
    end
end

% Set cathode parameters
for i=1:numel(cathode)
    S.([sideChar,'s1']).(['k',num2str(cathode(i))]).perc = 1/numel(cathode)*100;
    S.([sideChar,'s1']).(['k',num2str(cathode(i))]).pol = 1;
end

% Save stimulation parameters in template space output folder
save(fullfile(templateStimFolder, [patientName, '_desc-stimparameters.mat']), 'S');

% Save stimulation parameters in native space output folder
save(fullfile(nativeStimFolder, [patientName, '_desc-stimparameters.mat']), 'S');

% Fix missing atlasset in options
options.atlasset = 'DISTAL Nano (Ewert 2017)';

% Run VTA calculation
ea_genvat_horn([], S, sideInd, options, stimLabel);
