function json = ea_genrawimagesjson(BIDSRoot, subjId)
% [Re-]generate rawimages json file in case it's not present in subject's
% derivatives folder

warning('off', 'backtrace');
warning('Re-generating rawimages json file for "sub-%s" ...', subjId);

rawdataFolder = fullfile(BIDSRoot, 'rawdata', ['sub-', subjId]);

sessionFolders = flip(ea_regexpdir(rawdataFolder, '^ses-.*', 0, 'dir'));
sessions = regexp(sessionFolders, ['(?<=\', filesep, 'ses-)(.*)$'], 'match', 'once');

for s=1:length(sessions)
    datatypeFolders = ea_regexpdir(sessionFolders{s}, '[a-z]+', 0, 'dir');
    datatypes = regexp(datatypeFolders, ['(?<=\', filesep, ')([a-z]+)$'], 'match', 'once');
    for d=1:length(datatypes)
        niftiFiles = ea_regexpdir(datatypeFolders{d}, '.*\.nii(\.gz)?$', 0, 'file');
        for n=1:length(niftiFiles)
            parsed = parseBIDSFilePath(niftiFiles{n});
            if isfield(parsed, 'acq')
                key = [parsed.acq, '_', parsed.suffix];
            else
                key = parsed.suffix;
            end
            [~, json.(sessions{s}).(datatypes{d}).(key)] = ea_niifileparts(niftiFiles{n});
        end
    end
end

prefsFolder = fullfile(BIDSRoot, 'derivatives', 'leaddbs', ['sub-', subjId], 'prefs');
ea_mkdir(prefsFolder);
savejson('', json, fullfile(prefsFolder, ['sub-', subjId, '_desc-rawimages.json']));

warning('on', 'backtrace');
