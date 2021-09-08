function ea_coregpostopct(options)
% Entry function to coregister post-op CT

switch options.coregct.method
    case 'Advanced Normalization Tools (ANTs)'
        ea_coregctmri_ants(options);
    case 'Advanced Normalization Tools (ANTs), multiple runs'
        ea_coregctmri_ants_multiple(options);
    case 'Advanced Normalization Tools (ANTs), multiple runs + Subcortical Refine'
        ea_coregctmri_ants_multiple_refine(options);
    case 'Advanced Normalization Tools (ANTs) + Subcortical Refine'
        ea_coregctmri_ants_refine(options);
    case 'BRAINSFit'
        ea_coregctmri_brainsfit(options);
    case 'FSL FLIRT'
        ea_coregctmri_fsl(options);
end
