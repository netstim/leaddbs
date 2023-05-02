function ea_recentcallback(handles, type)
% Callback for recent patients/group analyses popupmenu

if get(handles.(['recent', type]), 'Value') == 1
    return
end

load([ea_getearoot, 'common', filesep, 'ea_recent', type, '.mat'],  'recentfolders');
if iscell(recentfolders)
    recentfolders = recentfolders(handles.(['recent', type]).Value-1);
end

if strcmp(['No recent ' type ' found'], recentfolders)
   return
end

if strcmp(type, 'datasets')
    datasetdir = handles.recentdatasets.String{handles.recentdatasets.Value};
    derivativesdir = fullfile(datasetdir, 'derivatives', 'leaddbs');
    rawdatadir = fullfile(datasetdir, 'rawdata');
    sourcedatadir = fullfile(datasetdir, 'sourcedata');
    if ~isfolder(datasetdir) || ~isfolder(derivativesdir) && ~isfolder(rawdatadir) && ~isfolder(sourcedatadir)
        ea_cprintf('CmdWinWarnin','Regenerating dataset folder structure since selected folder is missing or empty.\n');
        ea_mkdir(derivativesdir);
        ea_mkdir(rawdatadir);
        ea_mkdir(sourcedatadir);
    end
    ea_load_pts(handles, recentfolders);
elseif strcmp(type, 'patients')
    ea_load_pts(handles, recentfolders);
elseif strcmp(type, 'groups')
    ea_load_group(handles, recentfolders{1});
end
