function [file_path, status_path, releaseDir] = ea_input_programmer(options, numElectrodes)

    programmerDir = fullfile(options.earoot, 'programmer');

    file_path = fullfile(programmerDir, 'data.json');
    status_path = fullfile(programmerDir, 'status.json');
    input_file_path = fullfile(programmerDir, 'inputData.json');

    % Loop through each file path and create the file if it does not exist
    file_paths = {file_path, status_path, input_file_path};
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
    for i=1:numel(stimMatFile)
        stimJsonFile = strrep(stimMatFile{i}, '.mat', '.json');
        savejson('', load(stimMatFile{i}), stimJsonFile);
        fprintf('Converted and saved stimparameters.mat to %s\n', stimJsonFile);
    end

    inputStruct = struct;
    inputStruct.numElectrodes = numElectrodes;
    inputStruct.electrodeModel = options.elmodel;
    inputStruct.label = char(datetime('now', 'Format', 'yyyyMMddHHmmSS'));
    inputStruct.patientname = options.patientname;
    inputStruct.priorStims = {};
    if ~isempty(stimMatFile)
        inputStruct.priorStims = unique(fileparts(stimMatFile));
    end

    savejson('', inputStruct, input_file_path);

end