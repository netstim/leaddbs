function ea_initrecent(handles, type)
% Initialize recent patients/group analyses popupmenu

if ~exist('type', 'var')
    type = 'patients';
end

recentLog = [ea_getearoot, 'common', filesep, 'ea_recent', type, '.mat'];

if isfile(recentLog) && ismember('recentfolders', who('-file', recentLog))
    load(recentLog, 'recentfolders');
    recentfolders = recentfolders(isfolder(recentfolders));
else
    recentfolders = {['No recent ', type, ' found']};
end

save(recentLog, 'recentfolders');

% Refresh popupmenu
ea_refreshrecent(handles, type);
