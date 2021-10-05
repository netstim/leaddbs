function ea_refreshrecent(handles, type)
% Refresh recent patients/group analyses popupmenu

earoot = ea_getearoot;
load([earoot, 'common', filesep, 'ea_recent', type, '.mat'], 'recentfolders');

if ismember(['No recent ', type, ' found'], recentfolders) % No recent folder found, set to empty
    recentfolders = {};
else
    if strcmp(type, 'patients')
        recentfolders = regexp(recentfolders, ['(?<=\', filesep, 'leaddbs\', filesep, 'sub-).+'], 'match', 'once');
    elseif strcmp(type, 'groups')
        dataset = regexp(recentfolders, ['(?<=\', filesep, ')[^\', filesep, ']+', '(?=\', filesep, 'derivatives)'], 'match', 'once');
        analysis = regexp(recentfolders, ['(?<=\', filesep, 'leadgroup\', filesep, 'gs_).+'], 'match', 'once');
        recentfolders = strcat('dataset-', dataset, '_analysis-', analysis);
    end
end

recentfolders = [['Recent ', type, ':']; recentfolders];

% Refresh recent popupmenu, always set the value to 1
handles.recent.String = recentfolders;
handles.recent.Value = 1;
