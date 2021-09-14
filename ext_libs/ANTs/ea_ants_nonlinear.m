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
    disp(['Checking for slabs among structural images (assuming dominant structural file ',movingimage{end},' is a whole-brain acquisition)...']);
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
            sums(mov)=nan;
        end
    end
    if length(sums)>1 % multispectral warp
        is_slab = [sums(1:end-1) < (sums(end)*0.85) false];
    end
end

slab_present = any(is_slab);

if slab_present
    mnii.dt = [4,0];
    mnii.img = AllMX;
    tmaskdir = fullfile(fileparts(movingimage{1}), 'tmp');
    if ~exist(tmaskdir, 'dir')
        mkdir(tmaskdir);
    end
    mnii.fname = fullfile(tmaskdir, 'slabmask.nii');
    ea_write_nii(mnii);
    disp('Slabs found. Separating slabs to form an additional SyN stage.');
    slab_movingmask = mnii.fname;
else
    disp('No slabs found.');
    slab_movingmask = '';
end

% TODO: bids refactor
if options.prefs.machine.normsettings.ants_reinforcetargets
    if ~exist([options.root,options.patientname,filesep,'canat_combined.nii'],'file')
        ea_bet(movingimage{end},1,[options.root,options.patientname,filesep,'tmp.nii']);
        if ~exist([options.root,options.patientname,filesep,'tmp'],'dir')
            mkdir([options.root,options.patientname,filesep,'tmp']);
        end
        movefile([options.root,options.patientname,filesep,'tmp_mask.nii'],[options.root,options.patientname,filesep,'tmp',filesep,'brainmask.nii']);       
        ea_delete([options.root,options.patientname,filesep,'tmp.nii']);
        if slabspresent
            bmsk=ea_load_nii([options.root,options.patientname,filesep,'tmp',filesep,'brainmask.nii']);
            smsk=ea_load_nii([options.root,options.patientname,filesep,'tmp',filesep,'slabmask.nii']);
            smsk.img=smsk.img.*bmsk.img;
            smsk.fname=[options.root,options.patientname,filesep,'tmp',filesep,'brainslabmask.nii'];
            ea_write_nii(smsk);
            ea_combinenii_generic(movingimage,[options.root,options.patientname,filesep,'tmp',filesep,'brainslabmask.nii'],[options.root,options.patientname,filesep,'canat_combined.nii']);
        else
            ea_combinenii_generic(movingimage,[options.root,options.patientname,filesep,'tmp',filesep,'brainmask.nii'],[options.root,options.patientname,filesep,'canat_combined.nii']); % should get a brainmask here somehow.
        end
    end
end

% Load preset parameter set
apref = feval(eval(['@', options.prefs.machine.normsettings.ants_preset]), options.prefs.machine.normsettings);

% use fixed global correlations for fiducial helpers or segmentations
ccnsettg=options.prefs.machine.normsettings;
ccnsettg.ants_metric='Global Correlation';
ccpref = feval(eval(['@', options.prefs.machine.normsettings.ants_preset]), ccnsettg); 
ccpref.metric='MeanSquares';
ccpref.metricsuffix='';

metrics_prefix_suffix = cell(length(fixedimage),2);
[metrics_prefix_suffix{~is_segmentation,1}] = deal(apref.metric);
[metrics_prefix_suffix{~is_segmentation,2}] = deal(apref.metricsuffix);
[metrics_prefix_suffix{is_segmentation,1}] = deal(ccpref.metric);
[metrics_prefix_suffix{is_segmentation,2}] = deal(ccpref.metricsuffix);

for fi = 1:length(fixedimage)
    fixedimage{fi} = ea_path_helper(ea_niigz(fixedimage{fi}));
end
for fi = 1:length(movingimage)
    movingimage{fi} = ea_path_helper(ea_niigz(movingimage{fi}));
end

if length(fixedimage) ~= length(movingimage)
    ea_error('Please supply pairs of moving and fixed images (can be repetitive).');
end

outputimage = ea_path_helper(ea_niigz(outputimage));

basedir = [fileparts(mfilename('fullpath')), filesep];

if ispc
    HEADER = ea_path_helper([basedir, 'PrintHeader.exe']);
    ANTS = ea_path_helper([basedir, 'antsRegistration.exe']);
else
    HEADER = [basedir, 'PrintHeader.', computer('arch')];
    ANTS = [basedir, 'antsRegistration.', computer('arch')];
end

if ~ispc
    [~, imgsize] = system(['bash -c "', HEADER, ' ',fixedimage{1}, ' 2"']);
else
    [~, imgsize] = system([HEADER, ' ', fixedimage{1}, ' 2']);
end

imgsize = cellfun(@(x) str2double(x),ea_strsplit(imgsize,'x'));

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
    fixedmask=ea_path_helper([ea_space,'brainmask.nii.gz']);
else
    fixedmask='NULL';
end

% TODO: bids refactor
if false;exist(ea_niigz([outputdir,filesep,'mask_anatomy.nii']),'file')
    movingmask=ea_niigz([outputdir,filesep,'mask_anatomy.nii']);
else
    movingmask='NULL';
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
        ' --use-estimate-learning-rate-once ', ...
        ' --masks [',fixedmask,',',slab_movingmask,']'];

    for fi = 1:length(fixedimage)
        slabstage = [slabstage,...
            get_metric_command(fixedimage{fi}, movingimage{fi}, weights(fi), metrics_prefix_suffix(fi,:))];
    end
else
    slabstage = '';
end

% add subcortical refine stage:
if  options.prefs.machine.normsettings.ants_scrf
    synmaskconvergence = apref.convergence.scrf;
    synmaskshrinkfactors = apref.shrinkfactors.scrf;
    synmasksmoothingssigmas = apref.smoothingsigmas.scrf;

    synmaskstage = [' --transform ',apref.antsmode,apref.antsmode_suffix, ...
        ' --convergence ', synmaskconvergence, ...
        ' --shrink-factors ', synmaskshrinkfactors,  ...
        ' --smoothing-sigmas ', synmasksmoothingssigmas, ...
        ' --use-estimate-learning-rate-once ', ...
        ' --masks [',ea_path_helper([ea_space([],'subcortical'),'secondstepmask','.nii']),',',slab_movingmask,']'];
    for fi = 1:length(fixedimage)
        synmaskstage = [synmaskstage,...
            get_metric_command(fixedimage{fi}, movingimage{fi}, weights(fi), metrics_prefix_suffix(fi,:))];
    end
    
    % TODO: bids refactor
    strucs={'atlas'}; %{'STN','GPi','GPe','RN'};
    scnt=1;
    if false;options.prefs.machine.normsettings.ants_reinforcetargets
        for struc=1:length(strucs)
            disp(['Reinforcing ',strucs{scnt},' based on combined derived preop reconstruction']);
            synmaskstage = [synmaskstage,...
                ' --metric ',apref.metric,'[', ea_niigz([ea_space,strucs{struc}]), ',', ea_niigz([directory,'canat_combined']), ',',num2str(3),apref.metricsuffix,']'];
        end
    end
else
    synmaskstage = '';
end

ea_libs_helper
if options.prefs.machine.normsettings.ants_numcores
    setenv('ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS',options.prefs.machine.normsettings.ants_numcores) % no num2str needed since stored as string.
end

props.moving = movingimage{end};
props.fixed = fixedimage{end};
props.outputbase = outputbase;
props.ANTS = ANTS;
props.outputimage = outputimage;
props.rigidstage = rigidstage;
props.affinestage = affinestage;
props.synstage = synstage;
props.slabstage = slabstage;
props.synmaskstage = synmaskstage;
props.stagesep = options.prefs.machine.normsettings.ants_stagesep;

% if exist([fileparts(movingimage{1}),filesep,'glanatComposite.h5'],'file')
%     % clean old deformation field. this is important for cases where ANTs
%     % crashes and the user does not get an error back. Then, preexistant old transforms
%     % will be considered as new ones.
%     delete([fileparts(movingimage{1}),filesep,'glanatComposite.h5']);
% end

ea_submit_ants_nonlinear(props);

if exist('tmaskdir','var')
    ea_delete(tmaskdir);
end


function out = get_metric_command(fixed_image, moving_image, weight, metric_prefix_sufix)
    out = [' --metric ' metric_prefix_sufix{1} '[' fixed_image ',' moving_image ',' num2str(weight) metric_prefix_sufix{2} ']'];
 