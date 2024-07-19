function [file_path, releaseDir, status_path] = ea_input_programmer (options, M)

%     file_path = strcat(options.earoot, 'lead-dbs-programmer/data.json');
%     input_file_path = strcat(options.earoot, 'lead-dbs-programmer/inputData.json');

%   Initializing the communication files between matlab and electron
    file_path = strcat(options.earoot, 'programmergroup/data.json');
    status_path = strcat(options.earoot, 'programmergroup/status.json');
    input_file_path = strcat(options.earoot, 'programmergroup/inputData.json');
    fid = fopen(input_file_path, 'w');
    inputStruct = struct();

%     inputStruct.numElectrodes = length(elstruct.markers);

    inputStruct.electrodeModels = arrayfun(@(x) x.elmodel, M.elstruct, 'UniformOutput', false);
    inputStruct.label = M.guid;
    inputStruct.S = M.S;
    inputStruct.patientname = arrayfun(@(x) x.name, M.elstruct, 'UniformOutput', false);
    programmerDir = strcat(options.earoot, 'programmergroup');
    releaseDir = strcat(programmerDir, '/app/release/build');
%     stimDir = strcat(options.subj.stimDir, '/MNI152NLin2009bAsym');
%     stimFileName = strcat(options.patientname, '_desc-stimparameters.mat');
%     jsonFileName = strcat(options.patientname, '_desc-stimparameters.json');
%     directoryList = dir(stimDir);
%     % Initialize a cell array to store folder names
%     indicesToRemove = [];
%     for i = 1:numel(directoryList)
%         % Check if the item is a file, matches the name 'segmask.nii', or has the title 'Results_Rh'
%         if ~directoryList(i).isdir || strcmp(directoryList(i).name, 'segmask.nii') || strcmp(directoryList(i).name, 'Results_rh')
%             indicesToRemove(end+1) = i;
%         end
%     end
%     
%     % Remove the items at the specified indices
%     directoryList(indicesToRemove) = [];
%     if exist(stimDir, 'Dir')
%         inputStruct.priorStims = {};
%         inputStruct.priorStims = directoryList;
%     end
    fprintf(fid, '%s', jsonencode(inputStruct));
    fclose(fid);
%     for i=3:size(directoryList, 1)
%         stimLabel = directoryList(i).name;
%         currentStimDir = fullfile(stimDir, stimLabel);
%         matFileDir = fullfile(currentStimDir, stimFileName);
%         saveFileDir = fullfile(currentStimDir, jsonFileName);
%         % Load the MAT file
%         if exist(matFileDir, 'file')
%             matData = load(matFileDir);
%             
%             % Convert the MAT data to a JSON string
%             jsonData = jsonencode(matData);
%             
%             % Save the JSON string to a file
%             fid = fopen(saveFileDir, 'w');
%             if fid == -1
%                 error('Cannot open file for writing: %s', saveFileDir);
%             end
%             fprintf(fid, '%s', jsonData);
%             fclose(fid);
%             
%             fprintf('Converted and saved %s to %s\n', matFileDir, saveFileDir);
%         else
%             fprintf('MAT file does not exist: %s\n', matFileDir);
%         end
%     end

end