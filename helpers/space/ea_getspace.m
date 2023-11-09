function space=ea_getspace

prefsPath = ea_prefspath('mat');

if ~isfile(prefsPath)
    copyfile([ea_getearoot,'common',filesep,'ea_prefs_default.mat'], prefsPath, 'f');
end

load(prefsPath, 'machine');
space = machine.space;
