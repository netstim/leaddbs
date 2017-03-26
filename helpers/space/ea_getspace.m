function space=ea_getspace

%prefs=ea_prefs('');
home=ea_gethome;
if ~exist([home,'.ea_prefs.mat'],'file')
    copyfile([ea_getearoot,'common',filesep,'ea_prefs_default.mat'],[home,'.ea_prefs.mat']);
end
load([home,'.ea_prefs.mat']);
space=machine.space;