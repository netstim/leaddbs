function ea_restoreprefs(silent, cleanRecentHistory)

if ~exist('silent', 'var')
    silent = 0;
end

if ~exist('cleanRecentHistory', 'var')
    cleanRecentHistory = 0;
end

if ~silent
    answer = questdlg('Do you really want to reset your preferences to default?', 'Reset Lead-DBS preferences', 'Cancel', 'Sure', 'Cancel');
else
    answer = 'Sure';
end

if strcmp(answer,'Sure')
    copyfile([ea_getearoot, 'common', filesep, 'ea_prefs_default', ea_prefsext], ea_prefspath(ea_prefsext));
    copyfile([ea_getearoot, 'common', filesep, 'ea_prefs_default.mat'], ea_prefspath('mat'));
    ea_delete([ea_prefsdir, filesep, 'ea_ui.mat'])

    if cleanRecentHistory
        ea_delete([ea_prefsdir, filesep, 'ea_recentpatients.mat']);
        ea_delete([ea_prefsdir, filesep, 'ea_recentgroups.mat']);
    end
end
