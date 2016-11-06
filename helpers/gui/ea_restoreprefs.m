function ea_restoreprefs(varargin)

answer=questdlg('Do you really want to reset your preferences to default values?','Reset Lead-DBS preferences','Cancel','Sure','Cancel');
if strcmp(answer,'Sure')
    copyfile([ea_getearoot,'ea_prefs_default.m'],[ea_gethome,'.ea_prefs.m']);
end
