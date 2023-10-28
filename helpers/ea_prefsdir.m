function prefsDir = ea_prefsdir
% Return prefs directory

prefsDir = fullfile(ea_gethome, '.leaddbs');
ea_mkdir(prefsDir);
