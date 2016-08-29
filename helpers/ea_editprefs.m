function ea_editprefs(varargin)

if ~exist([ea_gethome,'.ea_prefs.m'],'file')
    copyfile([ea_getearoot,'ea_prefs_default.m'],[ea_gethome,'.ea_prefs.m']);
end
edit([ea_gethome,'.ea_prefs.m']);
