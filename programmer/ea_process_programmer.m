function [S] = ea_process_programmer(options)
    %% Initialize Paths, Read in data from programmer through data.json, Delete helper files (inputData, data)

    stimDir = fullfile(options.subj.stimDir, ea_getspace);

    importedS = loadjson(fullfile(stimDir, 'data.json'));
    
    ea_delete(fullfile(stimDir, 'data.json'));
    ea_delete(fullfile(stimDir, 'inputData.json'));

    if isfield(importedS, 'message')
        S = importedS;
        return;
    end

    %% Reorganize input for programmer to support lead-dbs convention

    S = importedS.S;
    amplitudeCell = arrayfun(@(i) S.amplitude(i, :), 1:size(S.amplitude, 1), 'UniformOutput', false);
    activeContactsCell = arrayfun(@(i) S.activecontacts(i, :), 1:size(S.activecontacts, 1), 'UniformOutput', false);
    S.amplitude = amplitudeCell;
    S.activecontacts = activeContactsCell;
    S.volume = [0 0];
    S.monopolarmodel = 0;
    S.sources=[1, 2, 3, 4];
end