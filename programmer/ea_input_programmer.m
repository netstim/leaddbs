function [file_path, releaseDir, input_file_path] = ea_input_programmer(options, numElectrodes)
% Prepare input for programmer

%% Handle output variables
programmerDir = fullfile(options.earoot, 'programmer');

file_path = fullfile(programmerDir, 'data.json');
if ~isfile(file_path)
    try
        savejson('', struct, file_path);
    catch ME
        error('Failed to create file: %s', file_path);
    end
end

releaseDir = fullfile(programmerDir, 'app', 'release');

%% Convert stimparameters.mat to json, handle inputData.json
inputStruct.patientname = options.patientname;
inputStruct.numElectrodes = numElectrodes;
inputStruct.electrodeModel = options.elmodel;
inputStruct.stimDir = fullfile(options.subj.stimDir, 'MNI152NLin2009bAsym');

stimDir = fullfile(options.subj.stimDir, ea_getspace);
ea_mkdir(stimDir);

stimFileName = [options.patientname, '_desc-stimparameters'];
stimMatFile = ea_regexpdir(stimDir, ['^', stimFileName, '\.mat$'], 1, 'f');

if isempty(stimMatFile)
    % Create new stimulation label, set S to empty
    inputStruct.labels = {[char(datetime('now', 'Format', 'yyyyMMddHHmmSS'))]};
    inputStruct.S = {};
else
    % Aggregate existing stimulations
    for i=1:numel(stimMatFile)
        load(stimMatFile{i}, 'S');

        if isempty(S.label)
            S.label = [char(datetime('now', 'Format', 'yyyyMMddHHmmSS'))];
        end
        inputStruct.labels{i} = S.label;

        if ~isfield(S, 'sources')
            S.sources = 1:4;
        end

        if ~isfield(S, 'volume')
            S.volume = [];
        end

        inputStruct.S(i) = S;

        stimJsonFile = strrep(stimMatFile{i}, '.mat', '.json');
        savejson('', S, stimJsonFile);
        fprintf('Converted and saved stimparameters.mat to %s\n', stimJsonFile);
    end
end

input_file_path = fullfile(stimDir, 'inputData.json');
savejson('', inputStruct, input_file_path);
