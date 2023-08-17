function ea_ants_run(cfg)
% Proxy to run ANTs registration based on provided configurations

% Make sure the transformation folder and output folder exist.
ea_mkdir(fileparts(cfg.outputbase));
ea_mkdir(fileparts(cfg.outputimage));

ants_transforms = dir(fullfile(fileparts(cfg.outputbase), '*_desc-ants.*'));

refinewarp = 0;
if ~isempty(ants_transforms) % prior ANTs transform found.
    if isfield(cfg, 'ants_usepreexisting')
        ants_usepreexisting = cfg.ants_usepreexisting;
    else
        prefs = ea_prefs;
        ants_usepreexisting = prefs.machine.normsettings.ants_usepreexisting;
    end
    switch ants_usepreexisting
        case 1 % ask
            answ = questdlg('We found existing ANTs transform files. Do you wish to build upon these transform (i.e. refine them) or discard them and start from scratch?','Old ANTs transform found.','Refine','Start from scratch','Start from scratch');
        case 2 % reuse
            answ = 'refine';
        case 3 % overwrite
            answ = 'start from scratch';
    end
    switch lower(answ)
        case 'refine'
            refinewarp = 1;
            cfg.rigidstage = '';
            cfg.affinestage = '';
        case 'start from scratch'
            % clean old deformation field. this is important for cases where ANTs
            % crashes and the user does not get an error back. Then, preexistant old transforms
            % will be considered as new ones.
            cellfun(@(x,y) ea_delete(fullfile(x,y)), {ants_transforms.folder}', {ants_transforms.name}')
            refinewarp = 0;
        otherwise
            return;
    end
end

fixedinit = fullfile(fileparts(fileparts(fileparts(cfg.moving))), 'masks', 'mask_template.nii');
if ~isfile(fixedinit)
    fixedinit = cfg.fixed;
end
fixedinit = ea_path_helper(fixedinit);

movinginit = fullfile(fileparts(fileparts(fileparts(cfg.moving))), 'masks', 'mask_anatomy.nii');
if ~isfile(movinginit)
    movinginit = cfg.moving;
end
movinginit = ea_path_helper(movinginit);

if refinewarp
    writecomposite = '0';
    forward_idx = cellfun(@(x) ~isempty(regexp(x, '.*from-anchorNative.*', 'once')), {ants_transforms.name});
    cfg.initial_transform = fullfile(ants_transforms(forward_idx).folder, ants_transforms(forward_idx).name);
    cfg.initial_inv_transform = fullfile(ants_transforms(~forward_idx).folder, ants_transforms(~forward_idx).name);
    initreg = [' --initial-moving-transform ', ea_path_helper(cfg.initial_transform)];
else
    writecomposite = '1';
    if isfield(cfg, 'initializationFeature') && ~isempty(cfg.initializationFeature)
        initializationFeature = cfg.initializationFeature;
    else
        % 0 for geometric center, 1 for image intensities, 2 for origin of the image
        initializationFeature = '0';
    end
    initreg = [' --initial-moving-transform [', fixedinit, ',', movinginit, ',', initializationFeature, ']'];
end

if isfield(cfg, 'histogrammatching') && ~isempty(cfg.histogrammatching)
    histogrammatching = cfg.histogrammatching;
else
    histogrammatching = '0';
end

if isfield(cfg, 'winsorize') && ~isempty(cfg.winsorize)
    winsorize = [' --winsorize-image-intensities [', cfg.winsorize, ']'];
else
    winsorize = '';
end

cmd = [cfg.ANTS, ' --verbose 1', ...
    ' --dimensionality 3', ...
    ' --float 1',...
    ' --write-composite-transform ', writecomposite, ...
    ' --output [',ea_path_helper(cfg.outputbase), ',', ea_path_helper(cfg.outputimage), ']', ...
    ' --interpolation Linear', ...
    ' --use-histogram-matching ', histogrammatching, ...
    winsorize, ...
    initreg, ...
    cfg.rigidstage, cfg.affinestage, cfg.synstage];

if isfield(cfg, 'slabstage')
    cmd = [cmd, cfg.slabstage];
end

if isfield(cfg, 'synmaskstage')
    cmd = [cmd, cfg.synmaskstage];
end

if isBIDSFileName(cfg.outputimage)
    logDir = fullfile(fileparts(fileparts(cfg.outputimage)), 'log');
    parsedStruct = parseBIDSFilePath(cfg.outputimage);
    ea_mkdir(logDir);
    antsCMDFile = [logDir, filesep, 'sub-', parsedStruct.sub, '_desc-antscmd.txt'];
else
    antsCMDFile = [fileparts(cfg.outputimage), filesep, 'ea_ants_command.txt'];
end

fid = fopen(antsCMDFile, 'a');
fprintf(fid, '%s:\n%s\n\n', char(datetime('now')), cmd);
fclose(fid);

status = ea_runcmd(cmd);

if status
    error(sprintf('ANTs normalization failed! Please check the log above for details.\nIn case it''s an out of memory error, reduce the number of threads in the ANTs settings might help.'));
end

if refinewarp
    ea_addrefinewarp(cfg);
else
    ea_conv_antswarps([cfg.outputbase, 'Composite.h5'], cfg.fixed, 'float');
    ea_conv_antswarps([cfg.outputbase, 'InverseComposite.h5'], cfg.moving, 'float');
end


function ea_addrefinewarp(props)

if exist([props.outputbase,'2Warp.nii.gz'],'file') % happens in second iteration of normalization refine
    refine_with_nth_prefix(props, '2');
elseif  exist([props.outputbase,'1Warp.nii.gz'],'file') % happens in second iteration of normalization refine
    refine_with_nth_prefix(props, '1');
end

% delete all old-version warps
ea_delete(props.initial_transform);
ea_delete(props.initial_inv_transform);
aux_warps = dir([props.outputbase '*Warp*']);
cellfun(@(x,y) ea_delete(fullfile(x,y)), {aux_warps.folder}', {aux_warps.name}');


function refine_with_nth_prefix(props, N)

outputformat = '.nii.gz';

applyTransforms = strrep(props.ANTS, 'antsRegistration', 'antsApplyTransforms');

cmd = [applyTransforms ' -r ' props.fixed ...
    ' -t '  ea_path_helper([props.outputbase N 'Warp.nii.gz']) ...
    ' -t '  ea_path_helper(props.initial_transform) ...
    ' -o [' ea_path_helper([props.outputbase 'Composite' outputformat]) ',1]' ...
    ' --float'];

invcmd = [applyTransforms ' -r ' props.moving ...
    ' -t '  ea_path_helper(props.initial_inv_transform) ...
    ' -t '  ea_path_helper([props.outputbase N 'InverseWarp.nii.gz']) ...
    ' -o [' ea_path_helper([props.outputbase 'InverseComposite' outputformat]) ',1]' ...
    ' --float'];

ea_runcmd(cmd);
ea_runcmd(invcmd);
