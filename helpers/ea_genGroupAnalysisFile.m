function analysisFile = ea_genGroupAnalysisFile(folder)
% Generate new group analysis file based on input dataset or group analysis folder.

if contains(folder, ['derivatives', filesep, 'leadgroup', filesep]) % Input is group analysis folder
    groupdir = fullfile(folder, filesep);
    [~, folderName] = fileparts(fileparts(groupdir));
    if ~startsWith(folderName, 'gs_')
        error('Group analysis folder should start with "gs_"!');
    end
    guid = erase(folderName, 'gs_');
elseif isfolder(fullfile(folder, 'derivatives')) % Input is dataset root folder
    groupFolderExists = 1; 
    while groupFolderExists
        guid = inputdlg('Set Group Analysis ID/Label [a-zA-Z0-9]:', '', 1, {datestr(datevec(now), 'yyyymmddHHMMSS')});
        if isempty(guid) % Cancel clicked, open inputdlg again
            continue;
        else
            guid = guid{1};
        end
        groupdir = fullfile(folder, 'derivatives', 'leadgroup', ['gs_', guid], filesep);
        groupFolderExists = isfolder(groupdir);
        if groupFolderExists % Folder already exist, need to specify a new one.
            waitfor(errordlg('Specified ID/Label already exists!', 'Error', 'modal'));
        end
    end
    mkdir(groupdir);
end

M = ea_initializeM(guid);
M.ui.groupdir = groupdir;

% Get dataset name
dataset = regexp(groupdir, ['(?<=\', filesep, ')[^\', filesep, ']+', '(?=\', filesep, 'derivatives)'], 'match', 'once');

analysisFile = [groupdir, 'dataset-', dataset, '_analysis-', guid, '.mat'];

save(analysisFile, 'M', '-v7.3');
