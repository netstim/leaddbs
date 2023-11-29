function ea_checkprefs
% Check prefs files and adapt for new prefs folder

prefsDir = ea_prefsdir;

prefs = ea_regexpdir(ea_gethome, '^\.ea_prefs.*', 0, 'f', 1);
if ~isempty(prefs)
    newprefs = replace(prefs, [ea_gethome, '.'], [prefsDir, filesep]);
    cellfun(@(x,y) movefile(x,y,'f'), prefs, newprefs);
end


prefs = ea_regexpdir(ea_prefsdir, '^ea_prefs\..*', 0, 'f', 1);
if ~isempty(prefs)
    newprefs = replace(prefs, 'ea_prefs.', 'ea_prefs_user.');
    cellfun(@(x,y) movefile(x,y,'f'), prefs, newprefs);
end

recent = ea_regexpdir([ea_getearoot, 'common'], '^\ea_recent.*', 0, 'f');
if ~isempty(recent)
    newrecent = replace(recent, [ea_getearoot, 'common'], prefsDir);
    cellfun(@(x,y) movefile(x,y,'f'), recent, newrecent);
end

ui = ea_regexpdir(ea_getearoot, '^\ea_ui.*', 0, 'f');
if ~isempty(ui)
    newui = replace(ui, ea_getearoot, [prefsDir, filesep]);
    cellfun(@(x,y) movefile(x,y,'f'), ui, newui);
end

extFolder = [ea_getearoot, 'ext_libs', filesep];

ea_delete({[extFolder, 'mambaforge']});

if isfolder(fullfile(ea_prefsdir, 'SlicerForLeadDBS'))
    ea_delete([extFolder, 'SlicerForLeadDBS']);
elseif isfolder([extFolder, 'SlicerForLeadDBS'])
    movefile([extFolder, 'SlicerForLeadDBS'], ea_prefsdir);
end

if isfolder(fullfile(ea_prefsdir, 'fastsurfer'))
    ea_delete([extFolder, 'fastsurfer', filesep, 'upstream']);
elseif isfolder([extFolder, 'fastsurfer', filesep, 'upstream'])
    movefile([extFolder, 'fastsurfer', filesep, 'upstream'], fullfile(ea_prefsdir, 'fastsurfer'));
end
