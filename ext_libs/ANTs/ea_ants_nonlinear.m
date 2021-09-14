function ea_ants_nonlinear(varargin)
% Wrapper for ANTs nonlinear registration

fixedimage = varargin{1};
movingimage = varargin{2};
outputimage = varargin{3};

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

if nargin >= 4
    weights = varargin{4};
else
    weights = ones(length(fixedimage),1);
end

if nargin >= 5
    options = varargin{5};
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

[outputdir, outputname, ~] = fileparts(outputimage);
if outputdir
    outputbase = [outputdir, filesep, outputname];
else
    outputbase = ['.', filesep, outputname];
end

if slabsupport
    disp(['Checking for slabs among structural images (assuming dominant structural file ',movingimage{end},' is a whole-brain acquisition)...']);
    for mov = 1:length(movingimage)
        if ~(weights(mov)>2.9) % segmentations
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
    slabspresent = 0; % default no slabs present.

    if length(sums)>1 % multispectral warp
        slabs = sums(1:end-1) < (sums(end)*0.85);
        if any(slabs) % one image is smaller than 0.7% of last (dominant) image, a slab is prevalent.
            slabmovingimage = ea_path_helper(movingimage(slabs)); % move slabs to new cell slabimage
            slabfixedimage = ea_path_helper(fixedimage(slabs));
            movingimage(slabs) = []; % remove slabs from movingimage
            fixedimage(slabs) = []; % remove slabs from movingimage

            % write out slab mask
            slabspresent = 1;
            mnii.dt = [4,0];
            mnii.img = AllMX;

            tmaskdir = fullfile(outputdir, 'tmp');
            if ~exist(tmaskdir, 'dir')
                mkdir(tmaskdir);
            end

            mnii.fname = [tmaskdir, filesep, 'slabmask.nii'];
            ea_write_nii(mnii);
            disp('Slabs found. Separating slabs to form an additional SyN stage.');
        else
            disp('No slabs found.');
        end
    end
else
    slabspresent = 0;
    impmasks = repmat({'nan'},length(movingimage),1);
end

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

directory = fileparts(movingimage{end});
if isempty(directory)
    directory = ['.', filesep];
else
    directory = [directory, filesep];
end

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
    if exist(ea_niigz([outputdir,filesep,'mask_template.nii']),'file')
        fixedmask=ea_niigz([outputdir,filesep,'mask_template.nii']);
    else
        fixedmask='NULL';
    end
end

if exist(ea_niigz([outputdir,filesep,'mask_anatomy.nii']),'file')
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
    if ~(weights(fi)>=3) % fiducial or segmentation, do not include here
        rigidstage = [rigidstage,...
            ' --metric ',apref.metric,'[', fixedimage{fi}, ',', movingimage{fi}, ',',num2str(weights(fi)),apref.metricsuffix,']'];
    end
end

affinestage = [' --transform Affine[0.15]'... % bit faster gradient step (see https://github.com/stnava/ANTs/wiki/Anatomy-of-an-antsRegistration-call)
    ' --convergence ', affineconvergence, ...
    ' --shrink-factors ', affineshrinkfactors ...
    ' --smoothing-sigmas ', affinesmoothingssigmas, ...
    ' --masks [',fixedmask,',',movingmask,']'];

for fi = 1:length(fixedimage)
    if ~(weights(fi)>=3) % fiducial or segmentation, do not include here
        affinestage = [affinestage,...
            ' --metric ',apref.metric,'[', fixedimage{fi}, ',', movingimage{fi}, ',',num2str(weights(fi)),apref.metricsuffix,']'];
    end
end

synstage = [' --transform ',apref.antsmode,apref.antsmode_suffix...
    ' --convergence ', synconvergence, ...
    ' --shrink-factors ', synshrinkfactors ...
    ' --smoothing-sigmas ', synsmoothingssigmas, ...
    ' --masks [',fixedmask,',',movingmask,']'];

for fi = 1:length(fixedimage)
    if weights(fi)>=3 % fiducial or segmentation
        if ~slabspresent % only add in fiducial / segmentation if no slabstage is added later.
            synstage = [synstage,...
                ' --metric ',ccpref.metric,'[', fixedimage{fi}, ',', movingimage{fi}, ',',num2str(weights(fi)),ccpref.metricsuffix,']'];
        end
    else
        synstage = [synstage,...
            ' --metric ',apref.metric,'[', fixedimage{fi}, ',', movingimage{fi}, ',',num2str(weights(fi)),apref.metricsuffix,']'];
    end
end

% Add slab stage
if slabspresent
    slabstage = [' --transform ',apref.antsmode,apref.antsmode_suffix...
        ' --convergence ', synconvergence, ...
        ' --shrink-factors ', synshrinkfactors ...
        ' --smoothing-sigmas ', synsmoothingssigmas, ...
        ' --use-estimate-learning-rate-once ', ...
        ' --masks [',fixedmask,',',ea_path_helper([tmaskdir,filesep,'slabmask.nii']),']'];
    fixedimage = [fixedimage,slabfixedimage];
    movingimage = [movingimage,slabmovingimage];

    for fi = 1:length(fixedimage)
        if weights(fi)>=3 % fiducial or segmentation
            slabstage = [slabstage,...
                ' --metric ',ccpref.metric,'[', fixedimage{fi}, ',', movingimage{fi}, ',',num2str(weights(fi)),ccpref.metricsuffix,']'];
        else
            slabstage = [slabstage,...
                ' --metric ',apref.metric,'[', fixedimage{fi}, ',', movingimage{fi}, ',',num2str(weights(fi)),apref.metricsuffix,']'];
        end
    end
else
    slabstage = '';
end

% add subcortical refine stage:
if  options.prefs.machine.normsettings.ants_scrf
    synmaskconvergence = apref.convergence.scrf;
    synmaskshrinkfactors = apref.shrinkfactors.scrf;
    synmasksmoothingssigmas = apref.smoothingsigmas.scrf;

    if slabspresent
        movingmask = ea_path_helper([tmaskdir,filesep,'slabmask.nii']);
    else
        movingmask = 'NULL';
    end

    synmaskstage = [' --transform ',apref.antsmode,apref.antsmode_suffix, ...
        ' --convergence ', synmaskconvergence, ...
        ' --shrink-factors ', synmaskshrinkfactors,  ...
        ' --smoothing-sigmas ', synmasksmoothingssigmas, ...
        ' --use-estimate-learning-rate-once ', ...
        ' --masks [',ea_path_helper([ea_space([],'subcortical'),'secondstepmask','.nii']),',',movingmask,']'];
    for fi = 1:length(fixedimage)
        if weights(fi)>=3 % fiducial or segmentation
            synmaskstage = [synmaskstage,...
                ' --metric ',ccpref.metric,'[', fixedimage{fi}, ',', movingimage{fi}, ',',num2str(weights(fi)),ccpref.metricsuffix,']'];
        else
            synmaskstage = [synmaskstage,...
                ' --metric ',apref.metric,'[', fixedimage{fi}, ',', movingimage{fi}, ',',num2str(weights(fi)),apref.metricsuffix,']'];
        end
    end
    
    strucs={'atlas'}; %{'STN','GPi','GPe','RN'};
    scnt=1;
    if options.prefs.machine.normsettings.ants_reinforcetargets
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
props.directory = directory;
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
