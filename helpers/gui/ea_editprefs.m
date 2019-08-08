function ea_editprefs(varargin)

if ~exist([ea_gethome,'.ea_prefs',ea_prefsext],'file')
    copyfile([ea_getearoot,'common',filesep,'ea_prefs_default',ea_prefsext],[ea_gethome,'.ea_prefs',ea_prefsext]);
end
if ~isdeployed
    edit([ea_gethome,'.ea_prefs.m']);
else
    try
        system(['subl ' ea_gethome,'.ea_prefs.json']); % TODO: add more programs to open file
    catch
        warning(['Unable to open. File location: ', ea_gethome, '.ea_prefs.json']);
    end
end