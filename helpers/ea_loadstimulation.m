function S = ea_loadstimulation(varargin)
% Load and check stimulation parameters

if nargin == 1
    stimFile = varargin{1};
else
    options = varargin{1};
    stimLabel = varargin{2};

    stimDir = fullfile(options.subj.stimDir, ea_nt(options), stimLabel);
    stimFile = fullfile(stimDir, ['sub-', options.subj.subjId, '_desc-stimparameters.mat']);
end

S = ea_checkStimParams(stimFile);
