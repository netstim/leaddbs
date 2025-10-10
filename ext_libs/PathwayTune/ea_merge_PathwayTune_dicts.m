function ea_merge_PathwayTune_dicts(stimFolder,mainDict,stimSetsDictPath_rh,stimSetsDictPath_lh,fixedSymtpomsDictPath,targetProfilesDictPath)
% Merge PathwayTune dictionaries into one master_dict.json
% By Butenko, konstantinmgtu@gmail.com

arguments
    stimFolder                  % path to the folder where master_dict.json to be stored 
    mainDict                    % path to netblend_dict_file
    stimSetsDictPath_rh         % path to StimSets_info of the right electrode, otherwise None
    stimSetsDictPath_lh         % path to StimSets_info of the left electrode, otherwise None
    fixedSymtpomsDictPath       % path to Fixed_symptoms
    targetProfilesDictPath      % path to target_profiles
end

%inputFilePath = '/home/forel/Documents/data/JB_project/fin_mSTN/netblend_dict_file.json';

% % 2. Read the entire JSON file content into a string
% fid = fopen(inputFilePath, 'r');
% if fid == -1
%     error('Could not open file %s for reading.', inputFilePath);
% end
% jsonStringFromFile = fread(fid, '*char')';
% fclose(fid);
% 
% % 3. Use jsondecode to convert the JSON string back into a MATLAB struct
% jsonDict = jsondecode(jsonStringFromFile);

if ~strcmp(stimSetsDictPath_rh,'None')
    % 2. Read the entire JSON file content into a string
    fid = fopen(stimSetsDictPath_rh, 'r');
    if fid == -1
        error('Could not open file %s for reading.', stimSetsDictPath_rh);
    end
    jsonStringFromFile = fread(fid, '*char')';
    fclose(fid);
    
    % 3. Use jsondecode to convert the JSON string back into a MATLAB struct
    StimSetsDict = jsondecode(jsonStringFromFile);
    mainDict.stim_sets_rh = StimSetsDict;
end

if ~strcmp(stimSetsDictPath_lh,'None')
    % 2. Read the entire JSON file content into a string
    fid = fopen(stimSetsDictPath_lh, 'r');
    if fid == -1
        error('Could not open file %s for reading.', stimSetsDictPath_lh);
    end
    jsonStringFromFile = fread(fid, '*char')';
    fclose(fid);
    
    % 3. Use jsondecode to convert the JSON string back into a MATLAB struct
    StimSetsDict = jsondecode(jsonStringFromFile);
    mainDict.stim_sets_lh = StimSetsDict;
end


fid = fopen(fixedSymtpomsDictPath, 'r');
if fid == -1
    error('Could not open file %s for reading.', fixedSymtpomsDictPath);
end
jsonStringFromFile = fread(fid, '*char')';
fclose(fid);
FixedSymptomsDict = jsondecode(jsonStringFromFile);
mainDict.fixed_symptom_weights = FixedSymptomsDict.fixed_symptom_weights;


fid = fopen(targetProfilesDictPath, 'r');
if fid == -1
    error('Could not open file %s for reading.', targetProfilesDictPath);
end
jsonStringFromFile = fread(fid, '*char')';
fclose(fid);
TargetProfilesDict = jsondecode(jsonStringFromFile);
mainDict.target_profiles = TargetProfilesDict;


jsonString = jsonencode(mainDict, 'PrettyPrint', true);
outputFilePath = [stimFolder,filesep,'master_dict.json'];
fid = fopen(outputFilePath, 'w');
if fid == -1
    error('Could not open file %s for writing.', outputFilePath);
end

% 4. Write the JSON string to the file
fprintf(fid, '%s', jsonString);

% 5. Close the file handle
fclose(fid);

disp(['Successfully saved formatted JSON to: ', outputFilePath]);
