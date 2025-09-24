function appFile = ea_getProgrammerPath
% Get the programmer path
% Install/update it when necessary

releaseDir = fullfile(options.earoot, 'programmer', 'app', 'release');
currentOS = ea_getarch;
zipFile = fullfile(releaseDir, ['LeadDBSProgrammer_', currentOS, '.zip']);

latestVersion = '1.1.0';
try
    installedVersion = loadjson(fullfile(ea_prefsdir, 'Programmer', 'version.json')).version;
catch
    installedVersion = '';
end

if ismac
    appFile = fullfile(ea_prefsdir, 'Programmer', 'LeadDBSProgrammer.app', 'Contents', 'MacOS', 'LeadDBSProgrammer');
    extractDir = fullfile(ea_prefsdir, 'Programmer');
elseif isunix
    appFile = fullfile(ea_prefsdir, 'Programmer', 'LeadDBSProgrammer', 'LeadDBSProgrammer');
    extractDir = fullfile(ea_prefsdir, 'Programmer', 'LeadDBSProgrammer');
else
    appFile = fullfile(ea_prefsdir, 'Programmer', 'LeadDBSProgrammer', 'LeadDBSProgrammer.exe');
    extractDir = fullfile(ea_prefsdir, 'Programmer', 'LeadDBSProgrammer');
end

if ~isfile(appFile) || ~strcmp(latestVersion, installedVersion)
    unzip(zipFile, extractDir);
    if ismac
        system(['xattr -cr ', ea_path_helper(fullfile(ea_prefsdir, 'Programmer', 'LeadDBSProgrammer.app'))]);
    end
    savejson('', latestVersion, fullfile(ea_prefsdir, 'Programmer', 'version.json'));
    savejson('', struct('LeadDBS_Path', ea_getearoot), fullfile(ea_prefsdir, 'Programmer', 'Preferences.json'));
end
