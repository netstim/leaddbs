function ea_coregpostopmr(options)
% Entry function to coregister post-op MRI to pre-op MRI

% Set fixed anchor image
anchor = options.subj.preopAnat.(options.subj.AnchorModality).coreg;

% Post-op MR modalities
postopModality = fieldnames(options.subj.postopAnat);

% Copy file tp preproc, take care of .nii.gz raw image
if strcmp(options.pref.niiFileExt, '.nii')
    cellfun(@(x) copyfile(options.subj.postopAnat.(x).raw, [options.subj.postopAnat.(x).preproc, '.gz']), postopModality, 'Uni', 0);
    cellfun(@(x) gunzip([options.subj.postopAnat.(x).preproc, '.gz']), postopModality, 'Uni', 0);
    cellfun(@(x) delete([options.subj.postopAnat.(x).preproc, '.gz']), postopModality, 'Uni', 0);
else
    cellfun(@(x) copyfile(options.subj.postopAnat.(x).raw, options.subj.postopAnat.(x).preproc), postopModality, 'Uni', 0);
end

% Set moving and output image
moving = cellfun(@(x) options.subj.postopAnat.(x).preproc, postopModality, 'Uni', 0);
output = cellfun(@(x) options.subj.postopAnat.(x).coreg, postopModality, 'Uni', 0);

% Check moving image existence
moving_exists = cellfun(@(x) isfile(x), moving);

% Check registration lock/approval status
output_approved = cellfun(@(x) logical(ea_reglocked(options, x)), output);

% Remove non-existing moving image and approved output image
moving(~moving_exists | output_approved) = [];
output(~moving_exists | output_approved) = [];

% Return if no image remains
if isempty(moving)
    return;
end

% Setup log
if options.prefs.diary
    ea_mkdir(fileparts(options.subj.coreg.log.logBaseName));
    diary([options.subj.coreg.log.logBaseName, 'MR', datestr(now, 'yyyymmddTHHMMss'), '.log']);
end

% Do coregistration
for i=1:length(moving)
    ea_dumpmethod(options, 'coreg', ea_getmodality(moving{i}));
    ea_coregimages(options, moving{i}, anchor, output{i}, [], options.prefs.mrcoreg.writeoutcoreg);
end

if options.prefs.diary
    diary off;
end
