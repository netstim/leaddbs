function [S] = ea_process_programmer(varargin)

if nargin == 1
    options = varargin{1};
    stimDir = fullfile(options.subj.stimDir, ea_getspace);
    try
        importedS = loadjson(fullfile(stimDir, 'data.json'));
    catch
        S = struct('message', 'Stimulation parameters not saved');
        return;
    end
    
    ea_delete(fullfile(stimDir, 'data.json'));
    ea_delete(fullfile(stimDir, 'inputData.json'));

    if isfield(importedS, 'message')
        S = importedS;
        return;
    end
elseif nargin == 2
    importedS = varargin{1};
end

S = importedS.S;
amplitudeCell = arrayfun(@(i) S.amplitude(i, :), 1:size(S.amplitude, 1), 'UniformOutput', false);
activeContactsCell = arrayfun(@(i) S.activecontacts(i, :), 1:size(S.activecontacts, 1), 'UniformOutput', false);
S.amplitude = amplitudeCell;
S.activecontacts = activeContactsCell;
S.volume = [0 0];
S.monopolarmodel = 0;
S.sources=[1, 2, 3, 4];

if ~isfield(S, 'template')
    S.template = 0;
end
