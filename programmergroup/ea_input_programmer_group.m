function [file_path, releaseDir, input_file_path] = ea_input_programmer_group (options, M)

%   Initializing the communication files between matlab and electron
    file_path = fullfile(M.root, 'data.json');
    input_file_path = fullfile(M.root, 'inputData.json');
    inputStruct = struct();
    inputStruct.electrodeModels = arrayfun(@(x) x.elmodel, M.elstruct, 'UniformOutput', false);
    inputStruct.label = strcat('gs_', M.guid);
    inputStruct.S = M.S;
    inputStruct.stimDir = M.root;
    inputStruct.patientname = arrayfun(@(x) x.name, M.elstruct, 'UniformOutput', false);
    programmerDir = strcat(options.earoot, 'programmergroup');
    releaseDir = strcat(programmerDir, '/app/release');
    savejson('', inputStruct, input_file_path);


end