function ea_ants_nonlinear_coreg(varargin)
% Wrapper for ANTs 3-stage nonlinear coregistration

fixedimage = varargin{1};
movingimage = varargin{2};
outputimage = varargin{3};

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

% Overwrite the default setting from GUI
normsettings.ants_preset = 'ea_antspreset_ants_default_synquick';

[outputdir, outputname, ~] = fileparts(outputimage);
if outputdir
    outputbase = [outputdir, filesep, outputname];
else
    outputbase = ['.', filesep, outputname];
end

% Load preset parameter set
apref = feval(eval(['@', normsettings.ants_preset]), normsettings);

directory = fileparts(movingimage);
if isempty(directory)
    directory = ['.', filesep];
else
    directory = [directory, filesep];
end

fixedimage = ea_path_helper(ea_niigz(fixedimage));
movingimage = ea_path_helper(ea_niigz(movingimage));

outputimage = ea_path_helper(ea_niigz(outputimage));

basedir = [fileparts(mfilename('fullpath')), filesep];

if ispc
    ANTS = ea_path_helper([basedir, 'antsRegistration.exe']);
else
    ANTS = [basedir, 'antsRegistration.', computer('arch')];
end

rigidstage = [' --transform Rigid[', apref.rigid.gradientstep, ']', ...
    ' --metric ', apref.rigid.metric, '[', fixedimage, ',', movingimage, ',', apref.rigid.metricparams, ']', ...
    ' --convergence ', apref.rigid.convergence, ...
    ' --shrink-factors ', apref.rigid.shrinkfactors, ...
    ' --smoothing-sigmas ', apref.rigid.smoothingsigmas, ...
    ' --masks [', fixedmask, ',', movingmask, ']'];

affinestage = [' --transform Affine[', apref.affine.gradientstep, ']', ...
    ' --metric ', apref.affine.metric, '[', fixedimage, ',', movingimage, ',', apref.affine.metricparams, ']', ...
    ' --convergence ', apref.affine.convergence, ...
    ' --shrink-factors ', apref.affine.shrinkfactors, ...
    ' --smoothing-sigmas ', apref.affine.smoothingsigmas, ...
    ' --masks [', fixedmask, ',', movingmask, ']'];

synstage = [' --transform SyN[', apref.syn.gradientstep, ']', ...
    ' --metric ', apref.syn.metric, '[', fixedimage, ',', movingimage, ',', apref.syn.metricparams, ']', ...
    ' --convergence ', apref.syn.convergence, ...
    ' --shrink-factors ', apref.syn.shrinkfactors, ...
    ' --smoothing-sigmas ', apref.syn.smoothingsigmas, ...
    ' --masks [', fixedmask, ',', movingmask, ']'];

ea_libs_helper;
if normsettings.ants_numcores
    setenv('ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS',normsettings.ants_numcores) % no num2str needed since stored as string.
end

props.ANTS = ANTS;
props.fixed = fixedimage;
props.moving = movingimage;
props.outputbase = outputbase;
props.outputimage = outputimage;
props.rigidstage = rigidstage;
props.affinestage = affinestage;
props.synstage = synstage;
props.directory = directory;
props.stagesep = normsettings.ants_stagesep;
props.ants_usepreexisting = normsettings.ants_usepreexisting; % Overwrite

ea_submit_ants_nonlinear(props);
