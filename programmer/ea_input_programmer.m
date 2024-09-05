function [file_path, releaseDir, input_file_path] = ea_input_programmer(options, numElectrodes)
% Prepare input for programmer

%% Handle output variables
programmerDir = fullfile(options.earoot, 'programmer');

file_path = fullfile(programmerDir, 'data.json');
input_file_path = fullfile(programmerDir, 'inputData.json');

% Loop through each file path and create the file if it does not exist
file_paths = {file_path, input_file_path};
for i = 1:length(file_paths)
    if ~isfile(file_paths{i})
        % Create an empty file
        fid = fopen(file_paths{i}, 'w');
        if fid == -1
            error('Could not create file: %s', file_paths{i});
        end
        fclose(fid);
    end
end

releaseDir = fullfile(programmerDir, 'app', 'release', 'build');

%% Convert stimparameters.mat to json, handle inputData.json
inputStruct.patientname = options.patientname;
inputStruct.numElectrodes = numElectrodes;
inputStruct.electrodeModel = options.elmodel;
inputStruct.stimDir = fullfile(options.subj.stimDir, 'MNI152NLin2009bAsym');

stimDir = fullfile(options.subj.stimDir, ea_getspace);
if ~exist(stimDir, 'dir')
    mkdir(stimDir); % Create the directory if it doesn't exist
end

stimFileName = [options.patientname, '_desc-stimparameters'];
stimMatFile = ea_regexpdir(stimDir, ['^', stimFileName, '\.mat$'], 1, 'f');

if isempty(stimMatFile)
    % Create new stimulation label, set S to empty
%     inputStruct.labels = {[char(datetime('now', 'Format', 'yyyyMMddHHmmSS')) ea_genid_rand(1,6)]};
    inputStruct.labels = {[char(datetime('now', 'Format', 'yyyyMMddHHmmSS'))]};
    inputStruct.S = {};
else
    % Aggregate existing stimulations
    for i=1:numel(stimMatFile)
        load(stimMatFile{i}, 'S');
        if isempty(S.label)
            % Generate a random label using datetime and a random ID
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
