function ea_addrecent(handles, uidir, type)
% Add new item to recent patients/group analyses

earoot = ea_getearoot;

load([earoot,'common',filesep,'ea_recent',type,'.mat'], 'recentfolders');
if strcmp(recentfolders, [ 'No recent ', type, ' found'])
    recentfolders = {};
end

recentfolders = unique([uidir'; recentfolders], 'stable');
if length(recentfolders) > 10
   recentfolders = recentfolders(1:10);
end

save([earoot, 'common', filesep, 'ea_recent', type, '.mat'], 'recentfolders');

try
    currentItem = recentfolders{handles.recent.Value};
catch
    currentItem = ['Recent ', type, ':'];
end
[~, idx] = ismember(currentItem, recentfolders);
ea_refreshrecent(handles, type, idx);
