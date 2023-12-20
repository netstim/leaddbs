function ea_switch2classic

LeadRoot = ea_getearoot;

if isfile([ea_prefsdir, filesep, 'ea_recentgroups.mat'])
    disp('Backup recent groups from develop branch ...');
    movefile([ea_prefsdir, filesep, 'ea_recentgroups.mat'], [ea_prefsdir, filesep, 'ea_recentgroups.mat.dev'])
end

if isfile([ea_prefsdir, filesep, 'ea_recentgroups.mat.classic'])
    disp('Restore recent groups from classic branch  ...');
    movefile([ea_prefsdir, filesep, 'ea_recentgroups.mat.classic'], [ea_prefsdir, filesep, 'ea_recentgroups.mat'])
end

if isfile([ea_prefsdir, filesep, 'ea_recentpatients.mat'])
    disp('Backup recent patients from develop branch ...');
    movefile([ea_prefsdir, filesep, 'ea_recentpatients.mat'], [ea_prefsdir, filesep, 'ea_recentpatients.mat.dev'])
end

if isfile([ea_prefsdir, filesep, 'ea_recentpatients.mat.classic'])
    disp('Restore recent patients from classic branch  ...');
    movefile([ea_prefsdir, filesep, 'ea_recentpatients.mat.classic'], [ea_prefsdir, filesep, 'ea_recentpatients.mat'])
end

if isfile([ea_prefsdir, filesep, 'ea_ui.mat'])
    disp('Backup ea_ui.mat from develop branch ...');
    movefile([ea_prefsdir, filesep, 'ea_ui.mat'], [ea_prefsdir, filesep, 'ea_ui.mat.dev'])
end

if isfile([ea_prefsdir, filesep, 'ea_ui.mat.classic'])
    disp('Restore ea_ui.mat from classic branch  ...');
    movefile([ea_prefsdir, filesep, 'ea_ui.mat.classic'], [ea_prefsdir, filesep, 'ea_ui.mat'])
end

if isfile(ea_prefspath)
    disp('Backup ea_prefs_user.m from develop branch ...');
    movefile(ea_prefspath,  ea_prefspath('.m.dev'))
end

if isfile(ea_prefspath('.m.classic'))
    disp('Restore ea_prefs_user.m from classic branch  ...');
    movefile(ea_prefspath('.m.classic'), ea_prefspath)
end

if isfile(ea_prefspath('mat'))
    disp('Backup ea_prefs_user.mat from develop branch ...');
    movefile(ea_prefspath('mat'),  ea_prefspath('.mat.dev'))
end

if isfile(ea_prefspath('.mat.classic'))
    disp('Restore ea_prefs_user.mat from classic branch  ...');
    movefile(ea_prefspath('.mat.classic'), ea_prefspath('mat'))
    % Fix possible wrong template string
    load(ea_prefspath('mat'), 'machine');
    machine.d2.backdrop = 'MNI_ICBM_2009b_NLIN_ASYM T1 (Fonov)';
    machine.togglestates.template = 'MNI_ICBM_2009b_NLIN_ASYM T1 (Fonov)';
    save(ea_prefspath('mat'), 'machine');
end

if isfile(ea_prefspath('json'))
    disp('Backup ea_prefs_user.json from develop branch  ...');
    movefile(ea_prefspath('json'),  ea_prefspath('.json.dev'))
end

if isfile(ea_prefspath('.json.classic'))
    disp('Restore ea_prefs_user.json from classic branch  ...');
    movefile(ea_prefspath('.json.classic'), ea_prefspath('json'))
end

disp('Switch LeadDBS branch to classic ...')
system(['git -C ', LeadRoot, ' stash']);
system(['git -C ', LeadRoot, ' checkout classic']);

ea_setpath;
rehash toolboxcache;
disp('LeadDBS search path updated.');
