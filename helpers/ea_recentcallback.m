function ea_recentcallback(handles, type)
% Callback for recent patients/group analyses popupmenu

if get(handles.recent, 'Value') == 1
    return
end

load([ea_getearoot, 'common', filesep, 'ea_recent', type, '.mat'],  'recentfolders');
if iscell(recentfolders)
    recentfolders = recentfolders(handles.recent.Value-1);
end

if strcmp(['No recent ' type ' found'], recentfolders)
   return
end

if strcmp(type, 'patients')
    ea_load_pts(handles, recentfolders);
elseif strcmp(type, 'groups')
    ea_load_group(handles, recentfolders{1});
end
