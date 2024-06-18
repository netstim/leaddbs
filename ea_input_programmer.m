function [file_path, releaseDir] = ea_input_programmer (options, elstruct)

    file_path = strcat(options.earoot, 'lead-dbs-programmer/data.json');
    input_file_path = strcat(options.earoot, 'lead-dbs-programmer/inputData.json');
    fid = fopen(input_file_path, 'w');
    inputStruct = struct();
    dt = datetime('now');
    % Convert the datetime object to a string in the format 'yyyymmddHHMMSS'
    formattedDate = datestr(dt, 'yyyymmddHHMMSS');
    inputStruct.numElectrodes = length(elstruct.markers);
    inputStruct.electrodeModel = options.elmodel;
    inputStruct.label = formattedDate;
    inputStruct.patientname = options.patientname;
    programmerDir = strcat(options.earoot, 'lead-dbs-programmer');
    releaseDir = strcat(programmerDir, '/release/build');
    stimDir = strcat(options.subj.stimDir, '/MNI152NLin2009bAsym');
    stimFileName = strcat(options.patientname, '_desc-stimparameters.mat');
    jsonFileName = strcat(options.patientname, '_desc-stimparameters.json');
    directoryList = dir(stimDir);
    % Initialize a cell array to store folder names
    indicesToRemove = [];
    for i = 1:numel(directoryList)
        % Check if the item is a file or matches the name 'segmask.nii'
        if ~directoryList(i).isdir || strcmp(directoryList(i).name, 'segmask.nii')
            indicesToRemove(end+1) = i;
        end
    end
    
    % Remove the items at the specified indices
    directoryList(indicesToRemove) = [];
    if exist(stimDir, 'Dir')
        inputStruct.priorStims = {};
        inputStruct.priorStims = directoryList;
    end
    fprintf(fid, '%s', jsonencode(inputStruct));
    fclose(fid);
end