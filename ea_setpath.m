function ea_setpath(options)
% Set or unset MATLAB search path for Lead-DBS
%
% When set the path, some subfolders (e.g. mambaforge) will be ignored.
% When unset the path, only Lead-DBS root folder will be kept, subfolders
% will be removed.

arguments
    options.unset {mustBeNumericOrLogical} = false
end

leaddbsRoot = ea_getearoot;
extDir = fullfile(leaddbsRoot, 'ext_libs');
fastsurferDir = fullfile(extDir, 'fastsurfer');
surficeDir = fullfile(extDir, 'surfice');

rootSubDirs = listDir(leaddbsRoot, {'ext_libs', 'release'});
extSubDirs = listDir(extDir, {'mambaforge', 'SlicerNetstim', 'SlicerForLeadDBS', '@'});
surficeSubDirs = listDir(surficeDir, {'surfice.app'});

dirToAdd = strjoin({leaddbsRoot, extDir, fastsurferDir, surficeDir}, pathsep);
subdirsToAdd = strjoin(cellfun(@genpath, [rootSubDirs; extSubDirs; surficeSubDirs], 'Uni', 0), pathsep);

finalPath = strjoin({dirToAdd, subdirsToAdd}, pathsep);

if options.unset
    finalPath = erase(finalPath, textBoundary('start') + leaddbsRoot);
    rmpath(finalPath);
else
    addpath(finalPath);
end

savepath;


function subDirs = listDir(inputDir, startsWithPattern)
% List subfolders and filter unwanted based on startsWithPattern
startsWithPattern = [startsWithPattern, '.'];

subDirs = dir(inputDir);
subDirs(~[subDirs.isdir]' | startsWith({subDirs.name}', startsWithPattern)) = [];
subDirs = fullfile({subDirs.folder}', {subDirs.name}');
