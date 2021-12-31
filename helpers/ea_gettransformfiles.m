function transformfiles=ea_gettransformfiles(options)

json = loadjson(options.subj.norm.log.method);

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
end

transformfiles.forward=[options.subj.norm.transform.forwardBaseName,warpSuffix];
transformfiles.inverse=[options.subj.norm.transform.inverseBaseName,warpSuffix];

