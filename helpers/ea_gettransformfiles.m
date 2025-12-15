function transformfiles=ea_gettransformfiles(options)

try
    json = loadjson(options.subj.norm.log.method);
catch
    % Lead-Connectome mode: Use SPM deformation field directly
    directory = [options.root, options.patientname, filesep];
    transformfiles.forward = fullfile(directory, 'y_ea_normparams.nii');
    transformfiles.inverse = fullfile(directory, 'y_ea_inv_normparams.nii');
    return;
end

if contains(json.method, 'ANTs')
    if isfield(json, 'custom') && json.custom
        % Custom full path of the transformation supplied.
        warpSuffix = '';
    elseif contains(json.method, 'affine')
        % Three-step affine normalization (Schonecker 2009) used
        warpSuffix = 'ants.mat';
    else
        warpSuffix = 'ants.nii.gz';
    end
elseif contains(json.method, 'FNIRT')
    warpSuffix='fnirt.nii.gz'; % correct?
elseif contains(json.method, 'SPM')
    warpSuffix='spm.nii';
elseif contains(json.method, 'EasyReg')
    warpSuffix='ants.nii.gz';
elseif contains(json.method, 'SynthMorph')
    warpSuffix='ants.nii.gz';
end

transformfiles.forward=[options.subj.norm.transform.forwardBaseName,warpSuffix];
transformfiles.inverse=[options.subj.norm.transform.inverseBaseName,warpSuffix];

