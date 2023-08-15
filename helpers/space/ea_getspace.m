function space=ea_getspace

home = ea_gethome;

if ~isfile([home,'.ea_prefs.mat'])
    copyfile([ea_getearoot,'common',filesep,'ea_prefs_default.mat'], [home,'.ea_prefs.mat'], 'f');
end

load([home, '.ea_prefs.mat'], 'machine');
space = machine.space;
