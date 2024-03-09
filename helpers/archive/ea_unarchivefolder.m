function ea_unarchivefolder(options)
% function to unarchive lead-dbs folders as lossless as possible.

% unzip sourcedata:
if exist([options.subj.sourcedataDir,'.zip'],'file')
    unzip([options.subj.sourcedataDir,'.zip'],fileparts(options.subj.sourcedataDir))
    ea_delete([options.subj.sourcedataDir,'.zip']);
end

% gunzip rawdata:
ea_recursiveniigunzip(options.subj.rawdataDir);
% gunzip preprocdir:
ea_recursiveniigunzip(options.subj.preprocDir);
% gunzip coregdir:
ea_recursiveniigunzip(options.subj.coregDir);
% gunzip brainshift:
ea_recursiveniigunzip(options.subj.brainshiftDir);
% gunzip reco:
ea_recursiveniigunzip(options.subj.reconDir);

% recover inverse transform and reapply normalization .nii
ea_recover_native_to_mni(options.subj.subjDir);
ea_normalize_apply_normalization(options);

function ea_recursiveniigunzip(folder)
flist=dir([folder,filesep,'**/*.gz']);
for f=1:length(flist)
gunzip(fullfile(flist(f).folder,flist(f).name));
ea_delete(fullfile(flist(f).folder,flist(f).name));
end



