function ea_setdefaults(itemname,value)

load([ea_getearoot,'common',filesep,'ea_prefs_default.mat']);
machine.(itemname)=value;
save([ea_getearoot,'common',filesep,'ea_prefs_default.mat'],'machine');
