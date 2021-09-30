function ea_refreshrecent(handles, type, idx)
% Refresh recent patients/group analyses popupmenu

earoot = ea_getearoot;
load([earoot, 'common', filesep, 'ea_recent', type, '.mat'], 'recentfolders');

if strcmp(type, 'patients')
    recentfolders = regexp(recentfolders, ['(?<=\', filesep, 'leaddbs\', filesep, 'sub-).+'], 'match', 'once');
elseif strcmp(type, 'groups')
    dataset = regexp(recentfolders, ['(?<=\', filesep, ')[^\', filesep, ']+', '(?=\', filesep, 'derivatives)'], 'match', 'once');
    analysis = regexp(recentfolders, ['(?<=\', filesep, 'leadgroup\', filesep, 'gs_).+'], 'match', 'once');
    recentfolders = strcat('dataset-', dataset, '_analysis-', analysis);
end

recentfolders = [['Recent ', type, ':']; recentfolders];

handles.recent.String = recentfolders;
if exist('idx', 'var')
   handles.recent.Value = idx+1;
end
