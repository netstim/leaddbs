function ea_submit_ants_nonlinear(props)

if props.stagesep
    ea_antsnl_multistep(props) % run each ANTs stage separately
else
    ea_antsnl_monostep(props) % run all ANTs stages together
end


function ea_antsnl_monostep(props)
directory=props.directory;
outputPrefix = strrep(props.outputbase, directory, '');
refinewarp=0;
if exist([props.outputbase,'Composite',ea_getantstransformext(directory)],'file') % prior ANTs transform found.
    if isfield(props, 'ants_usepreexisting')
        ants_usepreexisting = props.ants_usepreexisting;
    else
        prefs = ea_prefs;
        ants_usepreexisting = prefs.machine.normsettings.ants_usepreexisting;
    end

    switch ants_usepreexisting
        case 1 % ask
            answ=questdlg('We found existing ANTs transform files. Do you wish to build upon these transform (i.e. refine them) or discard them and start from scratch?','Old ANTs transform found.','Refine','Start from scratch','Start from scratch');
            switch lower(answ)
                case 'refine'
                    refinewarp=1;
                    props.rigidstage='';
                    props.affinestage='';
                case 'start from scratch'
                    ea_delete([props.outputbase,'Composite',ea_getantstransformext(directory)])
                    ea_delete([props.outputbase,'InverseComposite',ea_getantstransformext(directory)])
                    refinewarp=0;
            end
        case 2 % reuse
            refinewarp=1;
            props.rigidstage='';
            props.affinestage='';
        case 3 % overwrite
            % clean old deformation field. this is important for cases where ANTs
            % crashes and the user does not get an error back. Then, preexistant old transforms
            % will be considered as new ones.
            ea_delete([props.outputbase,'Composite',ea_getantstransformext(directory)])
            ea_delete([props.outputbase,'InverseComposite',ea_getantstransformext(directory)])
            refinewarp=0;
    end
end

if exist(ea_niigz([props.directory,filesep,'mask_template.nii']),'file')
    fixedinit=ea_niigz([props.directory,filesep,'mask_template.nii']);
else
    fixedinit=props.fixed;
end

if exist(ea_niigz([props.directory,filesep,'mask_anatomy.nii']),'file')
    movinginit=ea_niigz([props.directory,filesep,'mask_anatomy.nii']);
else
    movinginit=props.moving;
end

if refinewarp
    writecomposite = '0';
    initreg=[' --initial-moving-transform ',ea_path_helper([props.outputbase,'Composite',ea_getantstransformext(directory)])];
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

if isfield('props', 'slabstage')
    cmd = [cmd, props.slabstage];
end

if isfield('props', 'synmaskstage')
    cmd = [cmd, props.synmaskstage];
end

fid = fopen([props.directory,'ea_ants_command.txt'],'a');
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
    ea_addrefinewarp(props.directory);
end

reference = {props.fixed, props.moving};
ea_conv_antswarps(props.directory, outputPrefix, reference, '.nii.gz', 'float');


function ea_addrefinewarp(directory)

outputformat='.nii.gz';

basedir = [fileparts(mfilename('fullpath')), filesep];
if ispc
    applyTransforms = ea_path_helper([basedir, 'antsApplyTransforms.exe']);
else
    applyTransforms = [basedir, 'antsApplyTransforms.', computer('arch')];
end

template=ea_niigz([ea_space,'t1']);
[options.root,options.patientname]=fileparts(fileparts(directory)); % 'directory' is /a/b/c/
options.root=[options.root,filesep];
options.prefs=ea_prefs(options.patientname);
options=ea_assignpretra(options);
prenii=[directory,options.prefs.prenii_unnormalized];

if exist([directory,'glanat2Warp.nii.gz'],'file') % happens in second iteration of normalization refine
    cmd=[applyTransforms,' -r ',template,...
        ' -t ',ea_path_helper([directory,'glanat2Warp.nii.gz']),...
        ' -t ',ea_path_helper([directory,'glanatComposite',ea_getantstransformext(directory)]),...
        ' -o [',ea_path_helper([directory,'glanatComposite',outputformat]),',1]'];
    icmd=[applyTransforms,' -r ',ea_path_helper(prenii),...
        ' -t ',ea_path_helper([directory,'glanatInverseComposite',ea_getantstransformext(directory)]),...
        ' -t ',ea_path_helper([directory,'glanat2InverseWarp.nii.gz']),...
        ' -o [',ea_path_helper([directory,'glanatInverseComposite',outputformat]),',1]'];
elseif exist([directory,'glanat1Warp.nii.gz'],'file') % happens in third and upward iteration of normalization refine
    cmd=[applyTransforms,' -r ',template,...
        ' -t ',ea_path_helper([directory,'glanat1Warp.nii.gz']),...
        ' -t ',ea_path_helper([directory,'glanatComposite',ea_getantstransformext(directory)]),...
        ' -o [',ea_path_helper([directory,'glanatComposite',outputformat]),',1]'];
    icmd=[applyTransforms,' -r ',ea_path_helper(prenii),...
        ' -t ',ea_path_helper([directory,'glanatInverseComposite',ea_getantstransformext(directory)]),...
        ' -t ',ea_path_helper([directory,'glanat1InverseWarp.nii.gz']),...
        ' -o [',ea_path_helper([directory,'glanatInverseComposite',outputformat]),',1]'];
end

if exist('cmd','var')
    if ~ispc
        system(['bash -c "', cmd, '"']);
        system(['bash -c "', icmd, '"']);
    else
        system(cmd);
        system(icmd);
    end
end

% delete all old-version warps
ea_delete([directory,'glanat2InverseWarp.nii.gz']);
ea_delete([directory,'glanat1InverseWarp.nii.gz']);
ea_delete([directory,'glanat0Warp.nii.gz']);
ea_delete([directory,'glanat1Warp.nii.gz']);
ea_delete([directory,'glanat2Warp.nii.gz']);
ea_delete([directory,'glanat0GenericAffine.mat']);
ea_delete([directory,'glanatComposite.h5']);
ea_delete([directory,'glanatInverseComposite.h5']);


function ea_antsnl_multistep(props)
% Execute stages one for one

% linear:
cmd = [props.ANTS, ' --verbose 1', ...
    ' --dimensionality 3', ...
    ' --output ',ea_path_helper(props.outputbase),'', ...
    ' --interpolation Linear', ...
    ' --use-histogram-matching 0', ...
    ' --winsorize-image-intensities[0.005,0.995]', ...
    ' --float 1',...
    props.rigidstage, props.affinestage];

fid = fopen([props.directory,'ea_ants_command.txt'],'a');
fprintf(fid, '%s:\n%s\n\n', datestr(datetime('now')), cmd);
fclose(fid);

if ~ispc
    system(['bash -c "', cmd, '"']);
else
    system(cmd);
end

stack = [' --initial-moving-transform ',ea_path_helper(props.outputbase),'0GenericAffine.mat'];
tstack = {[' --transform [',ea_path_helper(props.outputbase),'0GenericAffine.mat,0]']};
itstack = {[' --transform [',ea_path_helper(props.outputbase),'0GenericAffine.mat,1]']};

% SyN:
cmd = [props.ANTS, ' --verbose 1', ...
    ' --dimensionality 3', ...
    stack ... % Apply Linear
    ' --output ',ea_path_helper(props.outputbase), 'Diff', ...
    ' --interpolation Linear', ...
    ' --use-histogram-matching 0', ...
    ' --winsorize-image-intensities[0.005,0.995]', ...
    ' --float 1',...
    props.synstage];

display(cmd)
fid = fopen([props.directory,'ea_ants_command.txt'],'a');
fprintf(fid, '%s:\n%s\n\n', datestr(datetime('now')), cmd);
fclose(fid);

if ~ispc
    system(['bash -c "', cmd, '"']);
else
    system(cmd);
end

stack = [' --initial-moving-transform ',ea_path_helper(props.outputbase),'Diff1Warp.nii.gz', ...
    ' --initial-moving-transform ',ea_path_helper(props.outputbase),'Diff0GenericAffine.mat'];
tstack = {[' --transform ',ea_path_helper(props.outputbase),'Diff1Warp.nii.gz'],...
    [' --transform [',ea_path_helper(props.outputbase),'Diff0GenericAffine.mat,0]']};
itstack = {[' --transform [',ea_path_helper(props.outputbase),'Diff0GenericAffine.mat,1]'],...
    [' --transform ',ea_path_helper(props.outputbase),'Diff1InverseWarp.nii.gz']};

% Slab:
if ~isempty(props.slabstage)
    cmd = [props.ANTS, ' --verbose 1', ...
        ' --dimensionality 3', ...
        stack ... % Apply Syn, Apply Linear
        ' --output ',ea_path_helper(props.outputbase), 'DiffSlab', ...
        ' --interpolation Linear', ...
        ' --use-histogram-matching 0', ...
        ' --winsorize-image-intensities[0.005,0.995]', ...
        ' --float 1',...
        props.slabstage];

    display(cmd)
    fid = fopen([props.directory,'ea_ants_command.txt'],'a');
    fprintf(fid, '%s:\n%s\n\n', datestr(datetime('now')), cmd);
    fclose(fid);

    if ~ispc
        system(['bash -c "', cmd, '"']);
    else
        system(cmd);
    end
    stack = [' --initial-moving-transform ',ea_path_helper(props.outputbase),'DiffSlab2Warp.nii.gz',...
        ' --initial-moving-transform ',ea_path_helper(props.outputbase),'DiffSlab1Warp.nii.gz',...
        ' --initial-moving-transform ',ea_path_helper(props.outputbase),'DiffSlab0GenericAffine.mat'];
    tstack = {[' --transform ',ea_path_helper(props.outputbase),'DiffSlab2Warp.nii.gz'],...
        [' --transform ',ea_path_helper(props.outputbase),'DiffSlab1Warp.nii.gz'],...
        [' --transform [',ea_path_helper(props.outputbase),'DiffSlab0GenericAffine.mat,0]']};
    itstack = {[' --transform [',ea_path_helper(props.outputbase),'DiffSlab0GenericAffine.mat,1]'],...
        [' --transform ',ea_path_helper(props.outputbase),'DiffSlab2InverseWarp.nii.gz']};
end

% SyNSubcorticalRefineStage:
if ~isempty(props.synmaskstage)
    cmd = [props.ANTS, ' --verbose 1', ...
        ' --dimensionality 3', ...
        stack ... % Apply Slab % Apply Syn % Apply Linear
        ' --output ',ea_path_helper(props.outputbase), 'DiffScrf', ...
        ' --interpolation Linear', ...
        ' --use-histogram-matching 0', ...
        ' --winsorize-image-intensities[0.005,0.995]', ...
        ' --float 1',...
        props.synmaskstage];

    display(cmd)
    fid = fopen([props.directory,'ea_ants_command.txt'],'a');
    fprintf(fid, '%s:\n%s\n\n', datestr(datetime('now')), cmd);
    fclose(fid);

    if ~ispc
        system(['bash -c "', cmd, '"']);
    else
        system(cmd);
    end
    stack = [' --initial-moving-transform ',ea_path_helper(props.outputbase),'DiffScrf2Warp.nii.gz',...
        ' --initial-moving-transform ',ea_path_helper(props.outputbase),'DiffScrf1Warp.nii.gz',...
        ' --initial-moving-transform ',ea_path_helper(props.outputbase),'DiffScrf0GenericAffine.mat'];
    tstack = {[' --transform ',ea_path_helper(props.outputbase),'DiffScrf2Warp.nii.gz'],...
        [' --transform ',ea_path_helper(props.outputbase),'DiffScrf1Warp.nii.gz'],...
        [' --transform [',ea_path_helper(props.outputbase),'DiffScrf0GenericAffine.mat,0]']};
    itstack = {[' --transform [',ea_path_helper(props.outputbase),'DiffScrf0GenericAffine.mat,1]'],...
        [' --transform ',ea_path_helper(props.outputbase),'DiffScrf2InverseWarp.nii.gz']};
end

%% combine all warps to a single one:
basedir = [fileparts(mfilename('fullpath')), filesep];
if ispc
    applyTransforms = ea_path_helper([basedir, 'antsApplyTransforms.exe']);
else
    applyTransforms = [basedir, 'antsApplyTransforms.', computer('arch')];
end

template = [ea_space,'t1.nii'];
outputformat = '.nii.gz';

tstring = [];
itstring = [];
for t = 1:length(tstack)-1
    tstring = [tstring,tstack{t}];
    itstring = [itstring,itstack{t}];
end
options = ea_getptopts(props.directory);
prenii = [props.directory,options.prefs.prenii_unnormalized];

% apply command:
cmd = [applyTransforms,' -r ',template,...
       tstring, ...
       ' -o [',ea_path_helper([props.directory,'glanatComposite',outputformat]),',1]'];
icmd = [applyTransforms,' -r ',ea_path_helper(prenii),...
        tstring, ...
        ' -o [',ea_path_helper([props.directory,'glanatInverseComposite',outputformat]),',1]'];

if exist('cmd','var')
    if ~ispc
        system(['bash -c "', cmd, '"']);
        system(['bash -c "', icmd, '"']);
    else
        system(cmd);
        system(icmd);
    end
end

% delete all old-version warps
warning('off')
ea_delete([props.directory,'glanatComposite.h5']);
ea_delete([props.directory,'glanatInverseComposite.h5']);
ea_delete([props.directory,'glanat0GenericeAffine.mat']);
ea_delete([props.directory,'glanat*Warp.mat']);
ea_delete([props.directory,'glanat*InverseWarp.mat']);
warning('on')
