function outputformat=ea_conv_antswarps(directory,outputformat)

if ~exist('outputformat','var')
    % check which version is present:
    if ~exist([directory,'glanatComposite.h5'],'file') && exist([directory,'glanatComposite.nii.gz'],'file')
        outputformat='.nii.gz';
    else
        outputformat='.h5';
    end
end

if ispc
    sufx='.exe';
else
    sufx=computer('arch');
end

antsApply=[ea_getearoot,'ext_libs',filesep,'ANTs',filesep,'antsApplyTransforms.',sufx];

template=[ea_space,'t2.nii'];
[options.root,options.patientname]=fileparts(fileparts(directory)); % 'directory' is /a/b/c/
options.root=[options.root,filesep];
options.prefs=ea_prefs(options.patientname);
options=ea_assignpretra(options);
prenii=[directory,options.prefs.prenii_unnormalized];

if exist([directory,'glanatComposite.h5'],'file') && ~strcmp(outputformat,'.h5')
    cmd=[antsApply,' -r ',template,...
        ' -t ',ea_path_helper([directory,'glanatComposite.h5']),...
        ' -o [',ea_path_helper([directory,'glanatComposite',outputformat]),',1]'];
    icmd=[antsApply,' -r ',ea_path_helper(prenii),...
        ' -t ',ea_path_helper([directory,'glanatInverseComposite.h5']),...
        ' -o [',ea_path_helper([directory,'glanatInverseComposite',outputformat]),',1]'];
elseif exist([directory,'glanatComposite.nii.gz'],'file') && ~strcmp(outputformat,'.nii.gz')
    cmd=[antsApply,' -r ',template,...
        ' -t ',ea_path_helper([directory,'glanatComposite.nii.gz']),...
        ' -o [',ea_path_helper([directory,'glanatComposite',outputformat]),',1]'];
    icmd=[antsApply,' -r ',ea_path_helper(prenii),...
        ' -t ',ea_path_helper([directory,'glanatInverseComposite.nii.gz']),...
        ' -o [',ea_path_helper([directory,'glanatInverseComposite',outputformat]),',1]'];
elseif exist([directory,'lanatComposite.h5'],'file')
    cmd=[antsApply,' -r ',template,...
        ' -t ',ea_path_helper([directory,'lanatComposite.h5']),...
        ' -o [',ea_path_helper([directory,'glanatComposite',outputformat]),',1]'];
    icmd=[antsApply,' -r ',ea_path_helper(prenii),...
        ' -t ',ea_path_helper([directory,'lanatInverseComposite.h5']),...
        ' -o [',ea_path_helper([directory,'glanatInverseComposite',outputformat]),',1]'];
elseif exist([directory,'glanat1Warp.nii.gz'],'file')
    cmd=[antsApply,' -r ',template,...
        ' -t [',ea_path_helper([directory,'glanat0GenericAffine.mat']),',0]',...
        ' -t ',ea_path_helper([directory,'glanat1Warp.nii.gz']),...
        ' -o [',ea_path_helper([directory,'glanatComposite',outputformat]),',1]'];
    icmd=[antsApply,' -r ',ea_path_helper(prenii),...
        ' -t ',ea_path_helper([directory,'glanat1InverseWarp.nii.gz']),...
        ' -t [',ea_path_helper([directory,'glanat0GenericAffine.mat']),',1]',...
        ' -o [',ea_path_helper([directory,'glanatInverseComposite',outputformat]),',1]'];
elseif exist([directory,'lanat1Warp.nii.gz'],'file')
    cmd=[antsApply,' -r ',template,...
        ' -t [',ea_path_helper([directory,'lanat0GenericAffine.mat']),',0]',...
        ' -t ',ea_path_helper([directory,'lanat1Warp.nii.gz']),...
        ' -o [',ea_path_helper([directory,'glanatComposite',outputformat]),',1]'];
    icmd=[antsApply,' -r ',ea_path_helper(prenii),...
        ' -t ',ea_path_helper([directory,'lanat1InverseWarp.nii.gz']),...
        ' -t [',ea_path_helper([directory,'lanat0GenericAffine.mat']),',1]',...
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
switch outputformat
    case '.nii.gz'
        ea_delete([directory,'glanatComposite.h5']);
        ea_delete([directory,'glanatInverseComposite.h5']);
    case '.h5'
        ea_delete([directory,'glanatComposite.nii.gz']);
        ea_delete([directory,'glanatInverseComposite.nii.gz']);
end
ea_delete([directory,'lanatComposite.h5']);
ea_delete([directory,'lanatInverseComposite.h5']);
ea_delete([directory,'glanat0GenericeAffine.mat']);
ea_delete([directory,'glanat1Warp.mat']);
ea_delete([directory,'glanat1InverseWarp.mat']);
ea_delete([directory,'lanat0GenericeAffine.mat']);
ea_delete([directory,'lanat1Warp.mat']);
ea_delete([directory,'lanat1InverseWarp.mat']);
