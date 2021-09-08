function ea_coregpreopmri(options)
% Entry function to coregister pre-op MRI
% ______________________________________________________________________________
% Copyright (C) 2015 Charite University Medicine Berlin, Movement Disorders Unit
% Andreas Horn

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
out_approved = cellfun(@(x) logical(ea_reglocked(options, x)), output);

% Remove non-existing moving image and approved output image
moving(~moving_exists | out_approved) = [];
output(~moving_exists | out_approved) = [];

% Return if no image remains
if isempty(moving)
    return;
end

% Setup log
ea_mkdir(fileparts(options.subj.coreg.log.logBaseName));
diary([options.subj.coreg.log.logBaseName, 'MR', datestr(now, 'yyyymmddTHHMMss'), '.log']);

for i=1:length(moving)
    ea_coregimages(options, moving{i}, anchor, output{i});

    % Reslice images if needed
    V1 = ea_open_vol(anchor);
    V2 = ea_open_vol(output{i});
    if ~isequal(V1.mat, V2.mat)
        ea_conformspaceto(anchor, output{i}, 1);
    end

    % Better slab support
    nii = ea_load_nii(V2.fname);
    nii.img(abs(nii.img)<0.0001) = 0;
    ea_write_nii(nii);
end

diary off;
