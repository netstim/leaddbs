function ea_generate_tpm

spmdir=[fileparts(which('spm')),filesep];
earoot=[fileparts(which('lead')),filesep];

switch spm('ver')
    case 'SPM8'
copyfile([spmdir,'toolbox',filesep,'Seg',filesep,'TPM.nii'],[earoot,'templates',filesep,'TPM.nii']);
    case 'SPM12'
copyfile([spmdir,'tpm',filesep,'TPM.nii'],[earoot,'templates',filesep,'TPM.nii']);
end

ea_reslice_nii([earoot,'templates',filesep,'TPM.nii'],[earoot,'templates',filesep,'TPM.nii'],[0.5 0.5 0.5],0,0,3);