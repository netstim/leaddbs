function ea_initrecent(handles, type)
% Initialize recent patients/group analyses popupmenu

if ~exist('type', 'var')
    type = 'patients';
end

earoot = ea_getearoot;

try
    load([earoot, 'common', filesep, 'ea_recent', type, '.mat'], 'recentfolders');
catch
    recentfolders = {['No recent ', type, ' found']};
end

save([earoot, 'common', filesep, 'ea_recent', type, '.mat'], 'recentfolders');

% Refresh popupmenu
ea_refreshrecent(handles, type);
