function space=ea_getspace

%prefs=ea_prefs('');
home = ea_gethome;
if ~isfile(ea_prefspath('mat'))
    copyfile([ea_getearoot,'common',filesep,'ea_prefs_default.mat'], ea_prefspath('mat'), 'f');
end
load(ea_prefspath('mat'));
space = machine.space;
