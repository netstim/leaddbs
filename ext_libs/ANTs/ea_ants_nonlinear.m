function ea_ants_nonlinear(fixedimage, movingimage, outputimage, weights, outputbase, options)
% Wrapper for ANTs nonlinear registration

if ischar(fixedimage)
    fixedimage = {fixedimage};
elseif ~iscell(fixedimage)
    ea_error('Please supply variable fixedimage as either char or cellstring');
end

if ischar(movingimage)
    movingimage = {movingimage};
elseif ~iscell(movingimage)
    ea_error('Please supply variable fixedimage as either char or cellstring');
end

if ~exist('weights', 'var') || isempty(weights)
    weights = ones(length(fixedimage),1);
end

if ~exist('outputbase', 'var') || isempty(outputbase)
    [directory, name, ~] = fileparts(outputimage);
    outputbase = fullfile(directory, name);
    clear directory name
end

if exist('options', 'var')
    options.prefs=ea_prefs; % refresh prefs in case called from recompute window with different settings.
else
    umachine = load([ea_gethome, '.ea_prefs.mat'], 'machine');
    options.prefs.machine.normsettings = umachine.machine.normsettings;
end

if isempty(which(options.prefs.machine.normsettings.ants_preset))
    fprintf(['\nCurrent ANTs preset "%s" is unavailable/deprecated!\n', ...
            'Will use the default preset "Effective: LowVariance, Default" instead.\n', ...
            'You may check your ANTs setting again later.\n\n'], ...
            options.prefs.machine.normsettings.ants_preset);
    % load default ANTs presets when current setting is unavailable/deprecated
    load([ea_getearoot, 'common', filesep, 'ea_prefs_default.mat'], 'machine');
    ants_default_preset = machine.normsettings.ants_preset;
    options.prefs.machine.normsettings.ants_preset = ants_default_preset;

    % save default ANTs presets to user preference file
    load([ea_gethome, '.ea_prefs.mat'], 'machine')
    machine.normsettings.ants_preset = ants_default_preset;
    save([ea_gethome, '.ea_prefs.mat'], 'machine', '-append');
end

slabsupport = 1; % check for slabs in anat files and treat slabs differently (add additional SyN stage only in which slabs are being used).

is_segmentation = weights >= 3;
is_slab = false(length(fixedimage),1);

if slabsupport
    if isBIDSFileName(movingimage{end})
        parsedStruct = parseBIDSFilePath(movingimage{end});
        anchorName = parsedStruct.suffix;
    else
        [~, anchorName] = ea_niifileparts(movingimage{end});
    end
    disp(['Checking for slabs among structural images (assuming anchor image ',anchorName,' is a whole-brain acquisition)...']);
    for mov = 1:length(movingimage)
        if ~is_segmentation(mov)
            mnii = ea_load_nii(movingimage{mov});
            mnii.img(abs(mnii.img)<0.0001)=nan;
            mnii.img=~isnan(mnii.img);
            if ~exist('AllMX','var')
                AllMX = mnii.img;
            else
                try
                    AllMX = AllMX.*mnii.img;
                catch
                    ea_error(sprintf('Multispectral acquisitions are not co-registered & resliced to anchor-modality. Please run co-registration first!\n%s', movingimage{mov}));
                end
            end
            sums(mov) = sum(mnii.img(:));
        else
            sums(mov) = nan;
        end
    end
    if length(sums)>1 % multispectral warp
        is_slab = [sums(1:end-1) < (sums(end)*0.85) false];
    end
end

slab_present = any(is_slab);

if slab_present
    mnii.dt(1) = 4;
    mnii.img = AllMX;
    tmaskdir = fullfile(fileparts(movingimage{1}), 'tmp');
    if ~exist(tmaskdir, 'dir')
        mkdir(tmaskdir);
    end
    mnii.fname = fullfile(tmaskdir, 'slabmask.nii');
    ea_write_nii(mnii);
    disp('Slabs found. Separating slabs to form an additional SyN stage.');
    slab_movingmask = ea_path_helper(mnii.fname);
else
    disp('No slabs found.');
    slab_movingmask = '';
end

% Load preset parameter set
apref = feval(eval(['@', options.prefs.machine.normsettings.ants_preset]), options.prefs.machine.normsettings);

% use fixed global correlations for fiducial helpers or segmentations
ccnsettg = options.prefs.machine.normsettings;
ccnsettg.ants_metric = 'Global Correlation';
ccpref = feval(eval(['@', options.prefs.machine.normsettings.ants_preset]), ccnsettg);
ccpref.metric = 'MeanSquares';
ccpref.metricsuffix = '';

metrics_prefix_suffix = cell(length(fixedimage),2);
[metrics_prefix_suffix{~is_segmentation,1}] = deal(apref.metric);
[metrics_prefix_suffix{~is_segmentation,2}] = deal(apref.metricsuffix);
[metrics_prefix_suffix{is_segmentation,1}] = deal(ccpref.metric);
[metrics_prefix_suffix{is_segmentation,2}] = deal(ccpref.metricsuffix);

for fi = 1:length(fixedimage)
    fixedimage{fi} = ea_niigz(fixedimage{fi});
end
for fi = 1:length(movingimage)
    movingimage{fi} = ea_niigz(movingimage{fi});
end

if length(fixedimage) ~= length(movingimage)
    ea_error('Please supply pairs of moving and fixed images (can be repetitive).');
end

outputimage = ea_niigz(outputimage);

basedir = [fileparts(mfilename('fullpath')), filesep];

ANTS = ea_getExec([basedir, 'antsRegistration'], escapePath = 1);


rigidconvergence = apref.convergence.rigid;
rigidshrinkfactors = apref.shrinkfactors.rigid;
rigidsmoothingssigmas = apref.smoothingsigmas.rigid;

affineconvergence = apref.convergence.affine;
affineshrinkfactors = apref.shrinkfactors.affine;
affinesmoothingssigmas = apref.smoothingsigmas.affine;

synconvergence = apref.convergence.syn;
synshrinkfactors = apref.shrinkfactors.syn;
synsmoothingssigmas = apref.smoothingsigmas.syn;

if options.prefs.machine.normsettings.ants_skullstripped
    fixedmask = ea_path_helper([ea_space,'brainmask.nii.gz']);
else
    fixedmask = 'NULL';
end

movingmask = fullfile(fileparts(fileparts(fileparts(movingimage{end}))), 'masks', 'mask_anatomy.nii');
if isfile(movingmask)
    movingmask = ea_path_helper(movingmask);
else
    movingmask = 'NULL';
end

rigidstage = [' --transform Rigid[0.25]' ... % bit faster gradient step (see https://github.com/stnava/ANTs/wiki/Anatomy-of-an-antsRegistration-call)
    ' --convergence ', rigidconvergence, ...
    ' --shrink-factors ', rigidshrinkfactors, ...
    ' --smoothing-sigmas ', rigidsmoothingssigmas, ...
    ' --masks [',fixedmask,',',movingmask,']'];

for fi = 1:length(fixedimage)
    if ~is_segmentation(fi) && ~is_slab(fi)
        rigidstage = [rigidstage,...
            get_metric_command(fixedimage{fi}, movingimage{fi}, weights(fi), metrics_prefix_suffix(fi,:))];
    end
end

affinestage = [' --transform Affine[0.15]'... % bit faster gradient step (see https://github.com/stnava/ANTs/wiki/Anatomy-of-an-antsRegistration-call)
    ' --convergence ', affineconvergence, ...
    ' --shrink-factors ', affineshrinkfactors ...
    ' --smoothing-sigmas ', affinesmoothingssigmas, ...
    ' --masks [',fixedmask,',',movingmask,']'];

for fi = 1:length(fixedimage)
    if ~is_segmentation(fi) && ~is_slab(fi)
        affinestage = [affinestage,...
            get_metric_command(fixedimage{fi}, movingimage{fi}, weights(fi), metrics_prefix_suffix(fi,:))];
    end
end

synstage = [' --transform ',apref.antsmode,apref.antsmode_suffix...
    ' --convergence ', synconvergence, ...
    ' --shrink-factors ', synshrinkfactors ...
    ' --smoothing-sigmas ', synsmoothingssigmas, ...
    ' --masks [',fixedmask,',',movingmask,']'];

for fi = 1:length(fixedimage)
    if is_segmentation(fi) && slab_present % slab stage with segmentation is added later
        continue
    else
        synstage = [synstage,...
            get_metric_command(fixedimage{fi}, movingimage{fi}, weights(fi), metrics_prefix_suffix(fi,:))];
    end
end

% Add slab stage
if slab_present
    slabstage = [' --transform ',apref.antsmode,apref.antsmode_suffix...
        ' --convergence ', synconvergence, ...
        ' --shrink-factors ', synshrinkfactors ...
        ' --smoothing-sigmas ', synsmoothingssigmas, ...
        ' --masks [',fixedmask,',',slab_movingmask,']'];

    for fi = 1:length(fixedimage)
        slabstage = [slabstage,...
            get_metric_command(fixedimage{fi}, movingimage{fi}, weights(fi), metrics_prefix_suffix(fi,:))];
    end
else
    slabstage = '';
end

% Add subcortical refine stage
if  options.prefs.machine.normsettings.ants_scrf
    synmaskconvergence = apref.convergence.scrf;
    synmaskshrinkfactors = apref.shrinkfactors.scrf;
    synmasksmoothingssigmas = apref.smoothingsigmas.scrf;

    synmaskstage = [' --transform ',apref.antsmode,apref.antsmode_suffix, ...
        ' --convergence ', synmaskconvergence, ...
        ' --shrink-factors ', synmaskshrinkfactors,  ...
        ' --smoothing-sigmas ', synmasksmoothingssigmas, ...
        ' --masks [',ea_path_helper([ea_space([],'subcortical'),'secondstepmask','.nii']),',',slab_movingmask,']'];
    for fi = 1:length(fixedimage)
        synmaskstage = [synmaskstage,...
            get_metric_command(fixedimage{fi}, movingimage{fi}, weights(fi), metrics_prefix_suffix(fi,:))];
    end
else
    synmaskstage = '';
end

ea_libs_helper
if options.prefs.machine.normsettings.ants_numcores
    setenv('ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS',options.prefs.machine.normsettings.ants_numcores)
end

cfg.moving = movingimage{end};
cfg.fixed = fixedimage{end};
cfg.outputbase = outputbase;
cfg.ANTS = ANTS;
cfg.outputimage = outputimage;
cfg.rigidstage = rigidstage;
cfg.affinestage = affinestage;
cfg.synstage = synstage;
cfg.slabstage = slabstage;
cfg.synmaskstage = synmaskstage;

ea_ants_run(cfg);

if exist('tmaskdir','var')
    ea_delete(tmaskdir);
end


function out = get_metric_command(fixed_image, moving_image, weight, metric_prefix_sufix)
    out = [' --metric ' metric_prefix_sufix{1} '[' ea_path_helper(fixed_image) ',' ea_path_helper(moving_image) ',' num2str(weight) metric_prefix_sufix{2} ']'];
