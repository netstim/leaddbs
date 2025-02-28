function analysisFile = ea_genGroupAnalysisFile(folder)
% Generate new group analysis file based on input dataset or group analysis folder.

folder = GetFullPath(folder);

if contains(folder, [filesep, 'derivatives', filesep, 'leadgroup', filesep]) % Input is group analysis folder
    groupdir = erase(folder, filesep + textBoundary("end"));
    guid = regexp(groupdir, ['(?<=\', filesep, 'leadgroup\', filesep, ').+'], 'match', 'once');
elseif isfolder(fullfile(folder, 'derivatives')) % Input is dataset root folder
    groupFolderExists = 1; 
    while groupFolderExists
        guid = inputdlg('Set Group Analysis ID/Label [a-zA-Z0-9]:', '', 1, {datestr(datevec(now), 'yyyymmddHHMMSS')});
        if isempty(guid) % Cancel clicked, open inputdlg again
            continue;
        else
            guid = guid{1};
        end
        groupdir = fullfile(folder, 'derivatives', 'leadgroup', guid);
        groupFolderExists = isfolder(groupdir);
        if groupFolderExists % Folder already exist, need to specify a new one.
            waitfor(errordlg('Specified ID/Label already exists!', 'Error', 'modal'));
        end
    end
end

ea_mkdir(groupdir);

M = ea_initializeM(guid);
M.root = fullfile(groupdir, filesep);

% Get dataset name
dataset = regexp(groupdir, ['(?<=\', filesep, ')[^\', filesep, ']+', '(?=\', filesep, 'derivatives)'], 'match', 'once');

analysisFile = fullfile(groupdir, ['dataset-', dataset, '_analysis-', guid, '.mat']);

save(analysisFile, 'M', '-v7.3');
