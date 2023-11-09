function ea_restoreprefs(varargin)

answer=questdlg('Do you really want to reset your preferences to default values?','Reset Lead-DBS preferences','Cancel','Sure','Cancel');
if strcmp(answer,'Sure')
    copyfile([ea_getearoot,'common',filesep,'ea_prefs_default',ea_prefsext],ea_prefspath(ea_prefsext));
    copyfile([ea_getearoot,'common',filesep,'ea_prefs_default.mat'],ea_prefspath('mat'));
end
