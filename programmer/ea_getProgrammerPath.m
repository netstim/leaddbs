function appFile = ea_getProgrammerPath
% Get the programmer path
% Install/update it when necessary

releaseDir = fullfile(ea_getearoot, 'programmer', 'app', 'release');
currentOS = ea_getarch;
zipFile = fullfile(releaseDir, ['LeadDBSProgrammer_', currentOS, '.zip']);

latestVersion = '0.1.2';
try
    installedVersion = loadjson(fullfile(ea_prefsdir, 'Programmer', 'version.json')).version;
catch
    installedVersion = '';
end

if ~isfile(zipFile) || ~strcmp(latestVersion, installedVersion)
    ea_delete(zipFile);
    downloadLatestRelease(['LeadDBSProgrammer_', currentOS, '.zip'], releaseDir)
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
    if ismac
        ea_delete(fullfile(extractDir, 'LeadDBSProgrammer.app'));
    else
        ea_delete(extractDir);
    end

    unzip(zipFile, extractDir);
    if ismac
        system(['xattr -cr ', ea_path_helper(fullfile(ea_prefsdir, 'Programmer', 'LeadDBSProgrammer.app'))]);
    end
    savejson('', struct('version', latestVersion), fullfile(ea_prefsdir, 'Programmer', 'version.json'));
    savejson('', struct('LeadDBS_Path', ea_getearoot), fullfile(ea_prefsdir, 'Programmer', 'Preferences.json'));
end


function downloadLatestRelease(assets, saveDir)

downloadURL = ['https://www.lead-dbs.org/release/download.php?id=programmer_' ea_getarch];

% Create the save directory if it doesn't exist
if ~isfolder(saveDir)
    mkdir(saveDir);
end

% Form the full destination path
destination = fullfile(saveDir, assets);

try
    % Download the file
    disp('Downloading LeadDBSProgrammer...');
    websave(destination, downloadURL, weboptions('Timeout', 10, 'ContentType', 'binary'));
    disp('Done.');
catch ME
    % Handle errors (e.g., 404, network issues)
    error('Failed to download "%s":\n%s', assets, ME.message);
end
