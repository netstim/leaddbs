function ea_addrecent(handles, uidir, type)
% Add new item to recent patients/group analyses

try
    load([ea_prefsdir, filesep, 'ea_recent',type,'.mat'], 'recentfolders');
catch
    recentfolders = {};
end

if ismember(['No recent ', type, ' found'], recentfolders)
    recentfolders = {};
end

if isrow(uidir)
    uidir = uidir';
end

recentfolders = unique(erase([uidir; recentfolders], filesep + textBoundary("end")), 'stable');
if length(recentfolders) > 10
   recentfolders = recentfolders(1:10);
end

save([ea_prefsdir, filesep, 'ea_recent', type, '.mat'], 'recentfolders');

ea_refreshrecent(handles, type);
