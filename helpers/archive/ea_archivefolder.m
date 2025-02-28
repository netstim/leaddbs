function ea_archivefolder(options)
% function to archive lead-dbs folders as lossless as possible to save
% storage space.

% zip sourcedata:
if ~exist([options.subj.sourcedataDir,'.zip'],'file')
    try zip([options.subj.sourcedataDir,'.zip'],options.subj.sourcedataDir); end
    ea_delete(options.subj.sourcedataDir);
end

% gz rawdata:
ea_recursiveniigz(options.subj.rawdataDir);
% gz preprocdir:
ea_recursiveniigz(options.subj.preprocDir);
% gz coregdir:
ea_recursiveniigz(options.subj.coregDir);
% gz brainshift:
ea_recursiveniigz(options.subj.brainshiftDir);
% gz reco:
ea_recursiveniigz(options.subj.reconDir);

% del normalization .nii
ea_delete(fullfile(options.subj.normDir,'anat','*.nii'));
ea_delete(fullfile(options.subj.normDir,'transformations',['sub-',options.subj.subjId,'_from-anchor*.nii.gz']));


function ea_recursiveniigz(folder)
flist=dir([folder,filesep,'**/*.nii']);
for f=1:length(flist)
gzip(fullfile(flist(f).folder,flist(f).name));
ea_delete(fullfile(flist(f).folder,flist(f).name));
end