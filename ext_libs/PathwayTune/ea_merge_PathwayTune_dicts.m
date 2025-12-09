function ea_merge_PathwayTune_dicts(stimFolder,mainDict,stimSetsDictPath_rh,stimSetsDictPath_lh,fixedSymtpomsDictPath,targetProfilesDictPath,masterDictPath)
% Merge PathwayTune dictionaries into one master_dict.json
% By Butenko, konstantinmgtu@gmail.com

arguments
    stimFolder                  % path to the folder where master_dict.json to be stored 
    mainDict                    % netblendict (pre-defined in ea_oss_optimizer.mlapp)
    stimSetsDictPath_rh         % path to StimSets_info of the right electrode, otherwise None
    stimSetsDictPath_lh         % path to StimSets_info of the left electrode, otherwise None
    fixedSymtpomsDictPath       % path to Fixed_symptoms
    targetProfilesDictPath      % path to target_profiles
    masterDictPath {mustBeText} = "None"
end

    % --- 1. Helper Function for Reading JSON ---
    function dataStruct = readJsonFile(filePath)
        % Reads a JSON file, decodes it, and returns the resulting structure.
        % Errors if the file cannot be read.
        if ~isfile(filePath)
            error('ea_merge:FileNotFound', 'JSON file not found: %s', filePath);
        end
        try
            jsonString = fileread(filePath);
            dataStruct = jsondecode(jsonString);
        catch ME
            error('ea_merge:JsonReadError', 'Error reading or decoding JSON from %s: %s', filePath, ME.message);
        end
    end

    % --- 3. Add StimSets for RH and LH Electrodes ---
    
    % Use `~strcmp(path, 'None')` or check `isfile` if "None" means the path is invalid/optional.
    % The original code used "None" as a sentinel value in string arguments.
    if ~strcmp(stimSetsDictPath_rh,'None')
        StimSetsDict_rh = readJsonFile(stimSetsDictPath_rh);
        mainDict.stim_sets_rh = StimSetsDict_rh;
    end

    if ~strcmp(stimSetsDictPath_lh,'None')
        StimSetsDict_lh = readJsonFile(stimSetsDictPath_lh);
        mainDict.stim_sets_lh = StimSetsDict_lh;
    end

    % --- 4. Handle Master Dictionary and Fixed Symptoms/Target Profiles ---
    
    masterDictExists = ~strcmp(masterDictPath,"None");

    if masterDictExists
        % Get fields from previously defined master_dict
        masterDict = readJsonFile(masterDictPath);
        
        % Always get target_profiles from masterDict if it exists
        if isfield(masterDict, 'target_profiles')
            mainDict.target_profiles = masterDict.target_profiles;
        end

        if isfile(fixedSymtpomsDictPath)
            % Use the new fixed symptom dict if available
            FixedSymptomsDict = readJsonFile(fixedSymtpomsDictPath);
            if isfield(FixedSymptomsDict, 'fixed_symptom_weights')
                mainDict.fixed_symptom_weights = FixedSymptomsDict.fixed_symptom_weights;
            end
        elseif isfield(masterDict, 'fixed_symptom_weights')
            % Otherwise, use the one from the previous master_dict
            mainDict.fixed_symptom_weights = masterDict.fixed_symptom_weights;
        end
    else
        % Master dict does not exist: Load fixed symptoms and target profiles from dedicated files
        
        % Load Fixed Symptoms
        FixedSymptomsDict = readJsonFile(fixedSymtpomsDictPath);
        if isfield(FixedSymptomsDict, 'fixed_symptom_weights')
            mainDict.fixed_symptom_weights = FixedSymptomsDict.fixed_symptom_weights;
        end
        
        % Load Target Profiles
        TargetProfilesDict = readJsonFile(targetProfilesDictPath);
        mainDict.target_profiles = TargetProfilesDict; % Assuming the whole dict is the target_profiles
    end

    % --- 5. Write the Master Dictionary to File ---
    outputFilePath = fullfile(stimFolder,'master_dict.json');

    % Use 'jsonencode' with 'PrettyPrint' for readability
    jsonString = jsonencode(mainDict, 'PrettyPrint', true);

    % Use 'filewrite' for simpler file output (available in recent MATLAB versions)
    % Alternatively, use fopen/fprintf/fclose with robust error handling (as below)
    
    fid = fopen(outputFilePath, 'w');
    if fid == -1
        error('ea_merge:FileWriteError', 'Could not open file %s for writing.', outputFilePath);
    end
    
    try
        fprintf(fid, '%s', jsonString);
        disp(['Successfully saved formatted JSON to: ', outputFilePath]);
    catch ME
        fclose(fid);
        rethrow(ME);
    end
    
    fclose(fid);

end