function ea_checkprefs
% Check prefs files and adapt for new prefs folder

prefsDir = ea_prefsdir;

prefs = ea_regexpdir(ea_gethome, '^\.ea_prefs.*', 0);
if ~isempty(prefs)
    newprefs = replace(prefs, [ea_gethome, '.'], [prefsDir, filesep]);
    cellfun(@(x,y) movefile(x,y,'f'), prefs, newprefs);
end

prefs = ea_regexpdir(ea_prefsdir, '^ea_prefs\..*', 0);
if ~isempty(prefs)
    newprefs = replace(prefs, 'ea_prefs.', 'ea_prefs_user.');
    cellfun(@(x,y) movefile(x,y,'f'), prefs, newprefs);
end

recent = ea_regexpdir([ea_getearoot, 'common'], '^\ea_recent.*', 0);
if ~isempty(recent)
    newrecent = replace(recent, [ea_getearoot, 'common'], prefsDir);
    cellfun(@(x,y) movefile(x,y,'f'), recent, newrecent);
end

ui = ea_regexpdir(ea_getearoot, '^\ea_ui.*', 0);
if ~isempty(ui)
    newui = replace(ui, ea_getearoot, [prefsDir, filesep]);
    cellfun(@(x,y) movefile(x,y,'f'), ui, newui);
end
