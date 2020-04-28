function ea_ants_nonlinear_coreg(varargin)
% Wrapper for ANTs 3-stage nonlinear coregistration

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

if nargin >= 4 && ~isempty(varargin{4})
    normsettings = varargin{4};
else
    umachine = load([ea_gethome, '.ea_prefs.mat'], 'machine');
    normsettings = umachine.machine.normsettings;
end

if nargin >= 5
    fixedmask = varargin{5};
else
    fixedmask = 'NULL';
end

if nargin >= 6
    movingmask = varargin{6};
else
    movingmask = 'NULL';
end

if isempty(which(normsettings.ants_preset))
    fprintf(['\nCurrent ANTs preset "%s" is unavailable/deprecated!\n', ...
            'Will use the default preset "Effective: LowVariance, Default" instead.\n', ...
            'You may check your ANTs setting again later.\n\n'], ...
            normsettings.ants_preset);
    % load default ANTs presets when current setting is unavailable/deprecated
    load([ea_getearoot, 'common', filesep, 'ea_prefs_default.mat'], 'machine');
    ants_default_preset = machine.normsettings.ants_preset;
    normsettings.ants_preset = ants_default_preset;

    % save default ANTs presets to user preference file
    load([ea_gethome, '.ea_prefs.mat'], 'machine')
    machine.normsettings.ants_preset = ants_default_preset;
    save([ea_gethome, '.ea_prefs.mat'], 'machine', '-append');
end

[outputdir, outputname, ~] = fileparts(outputimage);
if outputdir
    outputbase = [outputdir, filesep, outputname];
else
    outputbase = ['.', filesep, outputname];
end

% Load preset parameter set
apref = feval(eval(['@', normsettings.ants_preset]), normsettings);

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
    ANTS = ea_path_helper([basedir, 'antsRegistration.exe']);
else
    ANTS = [basedir, 'antsRegistration.', computer('arch')];
end

rigidconvergence = apref.convergence.rigid;
rigidshrinkfactors = apref.shrinkfactors.rigid;
rigidsmoothingssigmas = apref.smoothingsigmas.rigid;

affineconvergence = apref.convergence.affine;
affineshrinkfactors = apref.shrinkfactors.affine;
affinesmoothingssigmas = apref.smoothingsigmas.affine;

synconvergence = apref.convergence.syn;
synshrinkfactors = apref.shrinkfactors.syn;
synsmoothingssigmas = apref.smoothingsigmas.syn;

rigidstage = [' --transform Rigid[0.25]' ... % bit faster gradient step (see https://github.com/stnava/ANTs/wiki/Anatomy-of-an-antsRegistration-call)
    ' --convergence ', rigidconvergence, ...
    ' --shrink-factors ', rigidshrinkfactors, ...
    ' --smoothing-sigmas ', rigidsmoothingssigmas, ...
    ' --masks [',fixedmask,',',movingmask,']'];

for fi = 1:length(fixedimage)
    rigidstage = [rigidstage,...
        ' --metric ',apref.metric,'[', fixedimage{fi}, ',', movingimage{fi}, ',1',apref.metricsuffix,']'];
end

affinestage = [' --transform Affine[0.15]'... % bit faster gradient step (see https://github.com/stnava/ANTs/wiki/Anatomy-of-an-antsRegistration-call)
    ' --convergence ', affineconvergence, ...
    ' --shrink-factors ', affineshrinkfactors ...
    ' --smoothing-sigmas ', affinesmoothingssigmas, ...
    ' --masks [',fixedmask,',',movingmask,']'];

for fi = 1:length(fixedimage)
    affinestage = [affinestage,...
        ' --metric ',apref.metric,'[', fixedimage{fi}, ',', movingimage{fi}, ',1',apref.metricsuffix,']'];
end

synstage = [' --transform ',apref.antsmode,apref.antsmode_suffix...
    ' --convergence ', synconvergence, ...
    ' --shrink-factors ', synshrinkfactors ...
    ' --smoothing-sigmas ', synsmoothingssigmas, ...
    ' --masks [',fixedmask,',',movingmask,']'];


for fi = 1:length(fixedimage)
    synstage = [synstage,...
        ' --metric ',apref.metric,'[', fixedimage{fi}, ',', movingimage{fi}, ',1',apref.metricsuffix,']'];
end

ea_libs_helper;
if normsettings.ants_numcores
    setenv('ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS',normsettings.ants_numcores) % no num2str needed since stored as string.
end

props.moving = movingimage{end};
props.fixed = fixedimage{end};
props.outputbase = outputbase;
props.ANTS = ANTS;
props.outputimage = outputimage;
props.rigidstage = rigidstage;
props.affinestage = affinestage;
props.synstage = synstage;
props.directory = directory;
props.stagesep = normsettings.ants_stagesep;
props.ants_usepreexisting = normsettings.ants_usepreexisting; % Overwrite

ea_submit_ants_nonlinear(props);
