function ea_coregpreopmr(options)
% Entry function to coregister pre-op MRI

% Set anchor image
anchor = options.subj.coreg.anat.preop.(options.subj.AnchorModality);

% Set moving and output image
preopImage = rmfield(options.subj.preproc.anat.preop, options.subj.AnchorModality);
coregImage = rmfield(options.subj.coreg.anat.preop, options.subj.AnchorModality);
moving = struct2cell(preopImage);
output = struct2cell(coregImage);

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
ea_mkdir(fileparts(options.subj.coreg.log.logBaseName));
diary([options.subj.coreg.log.logBaseName, 'MR', datestr(now, 'yyyymmddTHHMMss'), '.log']);

% Do coregistration
for i=1:length(moving)
    ea_coregimages(options, moving{i}, anchor, output{i});

    % Better slab support
    nii = ea_load_nii(output{i});
    nii.img(abs(nii.img)<0.0001) = 0;
    ea_write_nii(nii);
end

ea_dumpmethod(options, 'coreg');

diary off;
