function ea_checkprefs
% Check prefs files and adapt new prefs folder

prefsDir = ea_prefsdir;

prefs = ea_regexpdir(ea_gethome, '^\.ea_prefs.*', 0, 'f', 1);

if ~isempty(prefs)
    newprefs = replace(prefs, [ea_gethome, '.'], [prefsDir, filesep]);
    cellfun(@(x,y) movefile(x,y), prefs, newprefs);
end
