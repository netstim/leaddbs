function done = ea_coregpreopmr(options)
% Entry function to coregister pre-op MRI

done = 0;

% Set anchor image
anchor = options.subj.coreg.anat.preop.(options.subj.AnchorModality);

% Set moving and output image
preopImage = rmfield(options.subj.preproc.anat.preop, options.subj.AnchorModality);
coregImage = rmfield(options.subj.coreg.anat.preop, options.subj.AnchorModality);

% FA is never coregistered on its own here: it is carried along with the B0
% (same ea_coregimages call, same transform), see below. Remove it from both
% structs to keep cell array sizes equal and avoid re-registration.
if isfield(preopImage, 'fa'), preopImage = rmfield(preopImage, 'fa'); end
if isfield(coregImage, 'fa'), coregImage = rmfield(coregImage, 'fa'); end
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
ea_mkdir(fullfile(options.subj.coregDir, 'transformations'));

for i=1:length(moving)
    ea_dumpmethod(options, 'coreg', ea_getmodality(moving{i}));

    % B0: apply the very same transform to the FA (DWI space) by passing it
    % as 'otherfiles'. This works for every method, nonlinear included.
    % The FA is processed in a temporary file and only renamed to its final
    % (anchor space) name once the B0 coregistration succeeded, so a failed
    % run never leaves a DWI-space FA under the anchor-space name.
    otherfiles = {};
    faCoreg = '';
    if strcmpi(ea_getmodality(moving{i}), 'b0')
        [faNative, faCoreg] = ea_getfapaths(options);
        if isfile(faNative)
            [faDir, faName, faExt] = fileparts(faCoreg);
            ea_mkdir(faDir);
            faTmp = fullfile(faDir, ['tmp_', faName, faExt]);
            ea_delete(faTmp);
            copyfile(faNative, faTmp);
            otherfiles = {faTmp};
        end
    end

    try
        affinefile = ea_coregimages(options, moving{i}, anchor, output{i}, otherfiles, 1);
    catch ME
        ea_delete(otherfiles);
        rethrow(ME);
    end

    if ~isempty(otherfiles)
        movefile(otherfiles{1}, faCoreg, 'f');
        otherfiles = {faCoreg};
        % Log FA like any coregistered modality: method.fa = method used
        % for the B0, approval.fa = 0 (the new FA has to be reviewed).
        ea_dumpmethod(options, 'coreg', 'fa');
        fprintf('FA coregistered with the B0->anchor transform: %s\n', faCoreg);
    end
    %% save Transforms
    switch lower(options.coregmr.method)
        case lower({'ANTs (Avants 2008)', 'ANTs'})
            movefile(affinefile{1},[options.subj.coreg.transform.(ea_getmodality(moving{i})).forwardBaseName, 'ants.mat']);
            movefile(affinefile{2},[options.subj.coreg.transform.(ea_getmodality(moving{i})).inverseBaseName, 'ants.mat']);
            % convert ANTS matrices to 4x4
            load([options.subj.coreg.transform.(ea_getmodality(moving{i})).forwardBaseName, 'ants.mat'])
            tmat = ea_antsmat2mat(AffineTransform_float_3_3,fixed);
            save([options.subj.coreg.transform.(ea_getmodality(moving{i})).forwardBaseName, 'ants44.mat'],'tmat')
            load([options.subj.coreg.transform.(ea_getmodality(moving{i})).inverseBaseName, 'ants.mat'])
            tmat = ea_antsmat2mat(AffineTransform_float_3_3,fixed);
            save([options.subj.coreg.transform.(ea_getmodality(moving{i})).inverseBaseName, 'ants44.mat'],'tmat')
        case lower({'FLIRT (Jenkinson 2001 & 2002)', 'FLIRT'})
            movefile(affinefile{1},[options.subj.coreg.transform.(ea_getmodality(moving{i})).forwardBaseName, 'flirt.mat']);
            movefile(affinefile{2},[options.subj.coreg.transform.(ea_getmodality(moving{i})).inverseBaseName, 'flirt.mat']);
            % convert affinefile from txt to tmat
            tmat = readmatrix([options.subj.coreg.transform.(ea_getmodality(moving{i})).forwardBaseName, 'flirt.mat'],'FileType','text');
            save([options.subj.coreg.transform.(ea_getmodality(moving{i})).forwardBaseName, 'flirt44.mat'],'tmat');
            tmat = readmatrix([options.subj.coreg.transform.(ea_getmodality(moving{i})).inverseBaseName, 'flirt.mat'],'FileType','text');
            save([options.subj.coreg.transform.(ea_getmodality(moving{i})).inverseBaseName, 'flirt44.mat'],'tmat');
        case lower({'SPM (Friston 2007)', 'SPM'})
            movefile(affinefile{1},[options.subj.coreg.transform.(ea_getmodality(moving{i})).forwardBaseName, 'spm.mat']);
            movefile(affinefile{2},[options.subj.coreg.transform.(ea_getmodality(moving{i})).inverseBaseName, 'spm.mat']);
            % also store tmat separatly analogous to the the other methods
            load([options.subj.coreg.transform.(ea_getmodality(moving{i})).forwardBaseName, 'spm.mat'],'tmat')
            save([options.subj.coreg.transform.(ea_getmodality(moving{i})).forwardBaseName, 'spm44.mat'],'tmat')
            load([options.subj.coreg.transform.(ea_getmodality(moving{i})).inverseBaseName, 'spm.mat'],'tmat')
            save([options.subj.coreg.transform.(ea_getmodality(moving{i})).inverseBaseName, 'spm44.mat'],'tmat')
        case lower({'BRAINSFit (Johnson 2007)', 'BRAINSFit'})
            movefile(affinefile{1},[options.subj.coreg.transform.(ea_getmodality(moving{i})).forwardBaseName, 'brainsfit.mat']);
            movefile(affinefile{2},[options.subj.coreg.transform.(ea_getmodality(moving{i})).inverseBaseName, 'brainsfit.mat']);
            % convert BRAINSFit matrices to 4x4
            load([options.subj.coreg.transform.(ea_getmodality(moving{i})).forwardBaseName, 'brainsfit.mat']);
            tmat = ea_antsmat2mat(AffineTransform_double_3_3,fixed);
            save([options.subj.coreg.transform.(ea_getmodality(moving{i})).forwardBaseName, 'brainsfit44.mat'],'tmat')
            load([options.subj.coreg.transform.(ea_getmodality(moving{i})).inverseBaseName, 'brainsfit.mat'])
            tmat = ea_antsmat2mat(AffineTransform_double_3_3,fixed);
            save([options.subj.coreg.transform.(ea_getmodality(moving{i})).inverseBaseName, 'brainsfit44.mat'],'tmat')
    end

    % Better slab support
    nii = ea_load_nii(output{i});
    nii.img(abs(nii.img)<0.0001) = 0;
    ea_write_nii(nii);

    % Write the FA after the B0, so the FA file is never older than the B0
    % it belongs to (ea_ensure_fa_and_fa2anat uses this to detect stale FA).
    for ofi = 1:length(otherfiles)
        if isfile(otherfiles{ofi})
            nii = ea_load_nii(otherfiles{ofi});
            nii.img(isnan(nii.img)) = 0; % resampling can leave NaNs (outside DWI field of view)
            nii.img(abs(nii.img)<0.0001) = 0;
            ea_write_nii(nii);
        end
    end
end

if options.prefs.diary
    diary off;
end

if options.overwriteapproved
    ea_segmask_cleanup(options);
    ea_cprintf('CmdWinWarnings', 'Preop MR coregistration has been rerun. Please also rerun normalization!\n');
end

done = 1;
