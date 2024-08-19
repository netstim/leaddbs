function [file_path, releaseDir] = ea_input_programmer(options, numElectrodes)

    programmerDir = fullfile(options.earoot, 'programmer');

    file_path = fullfile(programmerDir, 'data.json');
%     status_path = fullfile(programmerDir, 'status.json');
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
    
    stimDir = fullfile(options.subj.stimDir, ea_getspace);
    stimFileName = [options.patientname, '_desc-stimparameters'];

    stimMatFile = ea_regexpdir(stimDir, ['^', stimFileName, '\.mat$'], 1, 'f');
%     for i=1:numel(stimMatFile)
%         stimJsonFile = strrep(stimMatFile{i}, '.mat', '.json');
%         savejson('', load(stimMatFile{i}), stimJsonFile);
%         fprintf('Converted and saved stimparameters.mat to %s\n', stimJsonFile);
%     end

    labelsArray = {};
    exportStimFiles = struct();
%     exportStimFiles.S = {};
    for i = 1:numel(stimMatFile)
        % Load each .mat file
        load(stimMatFile{i});
        
        % Assuming S.label exists and is a cell array or a char array
        if isfield(S, 'label')
            labelsArray{i} = S.label; % Store labels in the array
        else
            labelsArray{i} = {}; % If no label, store an empty array
        end
        
        % Store the entire S struct in SStruct with a unique field name
        % Using the name of the file (without path and extension) as the field name
%         [~, name, ~] = fileparts(stimMatFile{i});
        exportStimFiles.S(i) = S;
        
        % Convert and save to a JSON file
        stimJsonFile = strrep(stimMatFile{i}, '.mat', '.json');
        savejson('', S, stimJsonFile);
        fprintf('Converted and saved %s to %s\n', stimMatFile{i}, stimJsonFile);
    end

    inputStruct = struct;
    inputStruct.numElectrodes = numElectrodes;
    inputStruct.electrodeModel = options.elmodel;
    inputStruct.label = char(datetime('now', 'Format', 'yyyyMMddHHmmSS'));
    if isempty(labelsArray)
        currentDate = datetime('now');
    
        % Extract year, month, and day
        year = currentDate.Year;
        month = sprintf('%02d', currentDate.Month); % Ensure two-digit month
        day = sprintf('%02d', currentDate.Day);     % Ensure two-digit day
        
        % Generate a random 6-digit number
        randomNums = randi([0 999999], 1);
        
        % Combine all parts into a string
        uniqueID = sprintf('%d%s%s%06d', year, month, day, randomNums);
        inputStruct.labels = {uniqueID};
    else
        
        inputStruct.labels = labelsArray;
    end
    inputStruct.patientname = options.patientname;
    inputStruct.stimDir = fullfile(options.subj.stimDir, 'MNI152NLin2009bAsym');
    if isfield(exportStimFiles, 'S')
        inputStruct.S = exportStimFiles.S;
    else
        inputStruct.S = {};
    end
%     if ~isempty(stimMatFile)
% %         inputStruct.priorStims = unique(fileparts(stimMatFile));
%         inputStruct.priorStims = fileparts(stimMatFile);
%     end

    savejson('', inputStruct, input_file_path);

end