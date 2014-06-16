%-----------------------------------------------------------------------
% Job configuration created by cfg_util (rev $Rev: 4252 $)
%-----------------------------------------------------------------------
matlabbatch{1}.spm.util.imcalc.input = {
                                        '/PA/Neuro/_projects/lead/lead/templates/segment/TPM.nii,1'
                                        '/PA/Neuro/_projects/lead/lead/templates/segment/STN_L_E.nii,1'
                                        };
matlabbatch{1}.spm.util.imcalc.output = 'tpm_STN_L_E.nii';
matlabbatch{1}.spm.util.imcalc.outdir = {''};
matlabbatch{1}.spm.util.imcalc.expression = 'i2';
matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
matlabbatch{1}.spm.util.imcalc.options.mask = 0;
matlabbatch{1}.spm.util.imcalc.options.interp = -4;
matlabbatch{1}.spm.util.imcalc.options.dtype = 4;


jobs{1}=matlabbatch;
cfg_util('run',jobs);
clear jobs matlabbatch