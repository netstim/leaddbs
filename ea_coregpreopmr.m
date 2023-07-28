function done = ea_coregpreopmr(options)
% Entry function to coregister pre-op MRI

done = 0;

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
if options.prefs.diary
    ea_mkdir(fileparts(options.subj.coreg.log.logBaseName));
    diary([options.subj.coreg.log.logBaseName, 'MR', datestr(now, 'yyyymmddTHHMMss'), '.log']);
end

% Do coregistration
for i=1:length(moving)
    ea_dumpmethod(options, 'coreg', ea_getmodality(moving{i}));

    affinefile = ea_coregimages(options, moving{i}, anchor, output{i},[],1);
    %% save Transforms
    switch lower(options.coregmr.method)
        case lower({'ANTs (Avants 2008)', 'ANTs'})
            movefile(affinefile{1},[options.subj.coreg.transform.(ea_getmodality(moving{i})).forwardBaseName, 'ants.mat']);
            movefile(affinefile{2},[options.subj.coreg.transform.(ea_getmodality(moving{i})).inverseBaseName, 'ants.mat']);
        case lower({'FLIRT (Jenkinson 2001 & 2002)', 'FLIRT'})
            movefile(affinefile{1},[options.subj.coreg.transform.(ea_getmodality(moving{i})).forwardBaseName, 'flirt.mat']);
            movefile(affinefile{2},[options.subj.coreg.transform.(ea_getmodality(moving{i})).inverseBaseName, 'flirt.mat']);
        case lower({'SPM (Friston 2007)', 'SPM'})
            movefile(affinefile{1},[options.subj.coreg.transform.(ea_getmodality(moving{i})).forwardBaseName, 'spm.mat']);
            movefile(affinefile{2},[options.subj.coreg.transform.(ea_getmodality(moving{i})).inverseBaseName, 'spm.mat']);
    end

    % Better slab support
    nii = ea_load_nii(output{i});
    nii.img(abs(nii.img)<0.0001) = 0;
    ea_write_nii(nii);
end

if options.prefs.diary
    diary off;
end

done = 1;
