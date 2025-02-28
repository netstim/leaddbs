function [input_file_path, releaseDir] = ea_input_programmer(options, numElectrodes)
% Prepare input for programmer

%% Handle output variables
releaseDir = fullfile(options.earoot, 'programmer', 'app', 'release');

%% Convert stimparameters.mat to json, handle inputData.json
inputStruct.nativeViewer = options.native;
inputStruct.estimateInTemplate = options.prefs.machine.vatsettings.estimateInTemplate;
inputStruct.patientname = options.patientname;
inputStruct.numElectrodes = numElectrodes;
inputStruct.electrodeModel = options.elmodel;

stimDir = fullfile(options.subj.stimDir, ea_getspace);
ea_mkdir(stimDir);

inputStruct.stimDir = stimDir;

stimFileName = [options.patientname, '_desc-stimparameters'];
stimMatFile = ea_regexpdir(stimDir, ['^', stimFileName, '\.mat$'], 1, 'f');
[~, stimFolder] = fileparts(fileparts(stimMatFile));
stimMatFile = stimMatFile(~startsWith(stimFolder, 'gs_'));

if isempty(stimMatFile)
    % Create new stimulation label, set S to empty
    inputStruct.labels = {[char(datetime('now', 'Format', 'yyyyMMddHHmmSS'))]};
    inputStruct.S(1) = ea_initializeS(inputStruct.labels, options);
else
    % Aggregate existing stimulations
    for i=1:numel(stimMatFile)
        S = ea_checkStimParams(stimMatFile{i});

        if isempty(S.label)
            S.label = [char(datetime('now', 'Format', 'yyyyMMddHHmmSS'))];
        end
        inputStruct.labels{i} = S.label;

        inputStruct.S(i) = S;

        stimJsonFile = strrep(stimMatFile{i}, '.mat', '.json');
        savejson('', S, stimJsonFile);
        fprintf('Converted and saved stimparameters.mat to %s\n', stimJsonFile);
    end
end

input_file_path = fullfile(stimDir, 'inputData.json');
savejson('', inputStruct, input_file_path);
