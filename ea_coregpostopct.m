function ea_coregpostopct(options)
% Entry function to coregister post-op CT to pre-op MRI

switch options.coregct.method
    case 'Advanced Normalization Tools (ANTs)'
        ea_coregpostopct_ants(options);
    case 'Advanced Normalization Tools (ANTs), multiple runs'
        ea_coregpostopct_ants_multiple(options);
    case 'Advanced Normalization Tools (ANTs), multiple runs + Subcortical Refine'
        ea_coregpostopct_ants_multiple_refine(options);
    case 'Advanced Normalization Tools (ANTs) + Subcortical Refine'
        ea_coregpostopct_ants_refine(options);
    case 'BRAINSFit'
        ea_coregpostopct_brainsfit(options);
    case 'FSL FLIRT'
        ea_coregpostopct_fsl(options);
end
