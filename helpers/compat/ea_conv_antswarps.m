function outputformat = ea_conv_antswarps(directory,prefix,reference,outputformat,float)

if ~exist('prefix','var') || isempty(prefix)
    prefix = 'glanat';
end

if ~exist('reference','var') || isempty(reference)
    fixedRef = [ea_space,'t1.nii'];
    [options.root,options.patientname] = fileparts(fileparts(directory)); % 'directory' is /a/b/c/
    options.root = [options.root,filesep];
    options.prefs = ea_prefs(options.patientname);
    options = ea_assignpretra(options);
    movingRef = [directory,options.prefs.prenii_unnormalized];
else
    fixedRef = reference{1};
    movingRef = reference{2};
end

if ~exist('outputformat','var')
    % check which version is present:
    if ~exist([directory,prefix,'Composite.h5'],'file') && exist([directory,prefix,'Composite.nii.gz'],'file')
        outputformat='.nii.gz';
    else
        outputformat='.h5';
    end
end

antsdir=[ea_getearoot,'ext_libs',filesep,'ANTs',filesep];
if ispc
    applyTransforms = ea_path_helper([antsdir, 'antsApplyTransforms.exe']);
else
    applyTransforms = [antsdir, 'antsApplyTransforms.', computer('arch')];
end

if exist([directory,prefix,'Composite.h5'],'file') && ~strcmp(outputformat,'.h5')
    cmd=[applyTransforms,' -r ',ea_path_helper(fixedRef),...
        ' -t ',ea_path_helper([directory,prefix,'Composite.h5']),...
        ' -o [',ea_path_helper([directory,prefix,'Composite',outputformat]),',1]'];
    icmd=[applyTransforms,' -r ',ea_path_helper(movingRef),...
        ' -t ',ea_path_helper([directory,prefix,'InverseComposite.h5']),...
        ' -o [',ea_path_helper([directory,prefix,'InverseComposite',outputformat]),',1]'];
elseif exist([directory,prefix,'Composite.nii.gz'],'file') && ~strcmp(outputformat,'.nii.gz')
    cmd=[applyTransforms,' -r ',ea_path_helper(fixedRef),...
        ' -t ',ea_path_helper([directory,prefix,'Composite.nii.gz']),...
        ' -o [',ea_path_helper([directory,prefix,'Composite',outputformat]),',1]'];
    icmd=[applyTransforms,' -r ',ea_path_helper(movingRef),...
        ' -t ',ea_path_helper([directory,prefix,'InverseComposite.nii.gz']),...
        ' -o [',ea_path_helper([directory,prefix,'InverseComposite',outputformat]),',1]'];
elseif exist([directory,prefix,'0GenericAffine.mat'],'file') && exist([directory,prefix,'1Warp.nii.gz'],'file')
    cmd=[applyTransforms,' -r ',ea_path_helper(fixedRef),...
        ' -t ',ea_path_helper([directory,prefix,'1Warp.nii.gz']),...
        ' -t [',ea_path_helper([directory,prefix,'0GenericAffine.mat']),',0]',...
        ' -o [',ea_path_helper([directory,prefix,'Composite',outputformat]),',1]'];
    icmd=[applyTransforms,' -r ',ea_path_helper(movingRef),...
        ' -t [',ea_path_helper([directory,prefix,'0GenericAffine.mat']),',1]',...
        ' -t ',ea_path_helper([directory,prefix,'1InverseWarp.nii.gz']),...
        ' -o [',ea_path_helper([directory,prefix,'InverseComposite',outputformat]),',1]'];
elseif exist([directory,'lanatComposite.h5'],'file')
    cmd=[applyTransforms,' -r ',ea_path_helper(fixedRef),...
        ' -t ',ea_path_helper([directory,'lanatComposite.h5']),...
        ' -o [',ea_path_helper([directory,prefix,'Composite',outputformat]),',1]'];
    icmd=[applyTransforms,' -r ',ea_path_helper(movingRef),...
        ' -t ',ea_path_helper([directory,'lanatInverseComposite.h5']),...
        ' -o [',ea_path_helper([directory,prefix,'InverseComposite',outputformat]),',1]'];
elseif exist([directory,'lanat1Warp.nii.gz'],'file')
    cmd=[applyTransforms,' -r ',ea_path_helper(fixedRef),...
        ' -t [',ea_path_helper([directory,'lanat0GenericAffine.mat']),',0]',...
        ' -t ',ea_path_helper([directory,'lanat1Warp.nii.gz']),...
        ' -o [',ea_path_helper([directory,prefix,'Composite',outputformat]),',1]'];
    icmd=[applyTransforms,' -r ',ea_path_helper(movingRef),...
        ' -t ',ea_path_helper([directory,'lanat1InverseWarp.nii.gz']),...
        ' -t [',ea_path_helper([directory,'lanat0GenericAffine.mat']),',1]',...
        ' -o [',ea_path_helper([directory,prefix,'InverseComposite',outputformat]),',1]'];
end

if exist('cmd','var')
    if exist('float', 'var')
        if ischar(float) && strcmp(float, 'float') || float
            cmd = [cmd, ' --float'];
            icmd = [icmd, ' --float'];
        end
    end

    cmd = [cmd, ' -v 1'];
    icmd = [icmd, ' -v 1'];

    if ~ispc
        system(['bash -c "', cmd, '"']);
        system(['bash -c "', icmd, '"']);
    else
        system(cmd);
        system(icmd);
    end
end

% delete all old-version warps
switch outputformat
    case '.nii.gz'
        ea_delete([directory,prefix,'Composite.h5']);
        ea_delete([directory,prefix,'InverseComposite.h5']);
    case '.h5'
        ea_delete([directory,prefix,'Composite.nii.gz']);
        ea_delete([directory,prefix,'InverseComposite.nii.gz']);
end

ea_delete([directory,prefix,'0GenericeAffine.mat']);
ea_delete([directory,prefix,'1Warp.mat']);
ea_delete([directory,prefix,'1InverseWarp.mat']);
ea_delete([directory,'lanatComposite.h5']);
ea_delete([directory,'lanatInverseComposite.h5']);
ea_delete([directory,'lanat0GenericeAffine.mat']);
ea_delete([directory,'lanat1Warp.mat']);
ea_delete([directory,'lanat1InverseWarp.mat']);
