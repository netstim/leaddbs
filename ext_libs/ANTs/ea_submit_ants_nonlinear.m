function ea_submit_ants_nonlinear(props)

if props.stagesep
    warning('ANTs multi step is depreciated, using mono step instead');
end

ea_antsnl_monostep(props) % run all ANTs stages together

function ea_antsnl_monostep(props)
ants_transforms = dir(fullfile(fileparts(props.outputbase), '*_desc-ants.*'));
refinewarp=0;
if ~isempty(ants_transforms) % prior ANTs transform found.
    if isfield(props, 'ants_usepreexisting')
        ants_usepreexisting = props.ants_usepreexisting;
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
            refinewarp=1;
            props.rigidstage='';
            props.affinestage='';
        case 'start from scratch'
            % clean old deformation field. this is important for cases where ANTs
            % crashes and the user does not get an error back. Then, preexistant old transforms
            % will be considered as new ones.
            cellfun(@(x,y) ea_delete(fullfile(x,y)), {ants_transforms.folder}', {ants_transforms.name}')
            refinewarp=0;
        otherwise
            return;
    end
end

% TODO: bids
if false;exist(ea_niigz([props.directory,filesep,'mask_template.nii']),'file')
    fixedinit=ea_niigz([props.directory,filesep,'mask_template.nii']);
else
    fixedinit=props.fixed;
end

% TODO: bids
if false;exist(ea_niigz([props.directory,filesep,'mask_anatomy.nii']),'file')
    movinginit=ea_niigz([props.directory,filesep,'mask_anatomy.nii']);
else
    movinginit=props.moving;
end

if refinewarp
    writecomposite = '0';
    forward_idx = cellfun(@(x) ~isempty(regexp(x, '.*from-anchorNative.*', 'once')), {ants_transforms.name});
    props.initial_transform = fullfile(ants_transforms(forward_idx).folder, ants_transforms(forward_idx).name);
    props.initial_inv_transform = fullfile(ants_transforms(~forward_idx).folder, ants_transforms(~forward_idx).name);
    initreg=[' --initial-moving-transform ', props.initial_transform];
else
    writecomposite = '1';
    if isfield(props, 'initializationFeature') && ~isempty(props.initializationFeature)
        initializationFeature = props.initializationFeature;
    else
         % 0 for geometric center, 1 for image intensities, 2 for origin of the image
        initializationFeature = '0';
    end
    initreg=[' --initial-moving-transform [', fixedinit, ',', movinginit, ',', initializationFeature, ']'];
end

if isfield(props, 'histogrammatching') && ~isempty(props.histogrammatching)
    histogrammatching = props.histogrammatching;
else
    histogrammatching = '0';
end

if isfield(props, 'winsorize') && ~isempty(props.winsorize)
    winsorize = [' --winsorize-image-intensities [', props.winsorize, ']'];
else
    winsorize = '';
end

cmd = [props.ANTS, ' --verbose 1', ...
    ' --dimensionality 3', ...
    ' --float 1',...
    ' --write-composite-transform ', writecomposite, ...
    ' --output [',ea_path_helper(props.outputbase), ',', props.outputimage, ']', ...
    ' --interpolation Linear', ...
    ' --use-histogram-matching ', histogrammatching, ...
    winsorize, ...
    initreg, ...
    props.rigidstage, props.affinestage, props.synstage];

if isfield(props, 'slabstage')
    cmd = [cmd, props.slabstage];
end

if isfield(props, 'synmaskstage')
    cmd = [cmd, props.synmaskstage];
end

% TODO: bidsify this name?
fid = fopen(fullfile(fileparts(fileparts(props.outputimage)), 'log', 'ea_ants_command.txt'), 'a');
fprintf(fid, '%s:\n%s\n\n', datestr(datetime('now')), cmd);
fclose(fid);

if ~ispc
    status=system(['bash -c "', cmd, '"']);
else
    status=system(cmd);
end

if status
   ea_error('ANTs normalization failed - likely due to out of memory problems. Please try a different normalization strategy or reduce the number of threads in the ANTs settings dialogue.');
end

if refinewarp
    ea_addrefinewarp(props);
else
    ea_conv_antswarps([props.outputbase, 'Composite.h5'], props.fixed, 'float');
    ea_conv_antswarps([props.outputbase, 'InverseComposite.h5'], props.moving, 'float');
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

outputformat='.nii.gz';

applyTransforms = strrep(props.ANTS, 'antsRegistration', 'antsApplyTransforms');

cmd = [applyTransforms ' -r ' props.fixed ...
    ' -t '  ea_path_helper([props.outputbase N 'Warp.nii.gz']) ...
    ' -t '  ea_path_helper(props.initial_transform) ...
    ' -o [' ea_path_helper([props.outputbase 'Composite' outputformat]) ',1]' ...
    ' --float'];

icmd = [applyTransforms ' -r ' props.moving ...
    ' -t '  ea_path_helper(props.initial_inv_transform) ...
    ' -t '  ea_path_helper([props.outputbase N 'InverseWarp.nii.gz']) ...
    ' -o [' ea_path_helper([props.outputbase 'InverseComposite' outputformat]) ',1]' ...
    ' --float'];

if ~ispc
    system(['bash -c "', cmd, '"']);
    system(['bash -c "', icmd, '"']);
else
    system(cmd);
    system(icmd);
end