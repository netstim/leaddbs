function ea_apply_normalization(options)
% Wrapper to apply normalization to the unnormalized files.
% __________________________________________________________________________________
% Copyright (C) 2014 Charite University Medicine Berlin, Movement Disorders Unit
% Andreas Horn

json = loadjson(options.subj.norm.log.method);

if contains(json.method, {'ANTs', 'EasyReg', 'SynthMorph', 'SPM'})
    ea_ants_apply_transforms(options);
elseif contains(json.method, 'FNIRT')
    ea_fsl_apply_normalization(options);
elseif contains(json.method, 'SPM')
    % Convert SPM deformation field to ITK format when necessary
    ea_convert_spm_warps(options.subj);
    ea_ants_apply_transforms(options);
end
