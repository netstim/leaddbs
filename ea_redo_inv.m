function ea_redo_inv(directory,options)
matlabbatch{1}.spm.util.defs.comp{1}.inv.comp{1}.def = {[directory,'y_ea_normparams.nii']};
matlabbatch{1}.spm.util.defs.comp{1}.inv.space = {[directory,options.prefs.prenii_unnormalized]};
matlabbatch{1}.spm.util.defs.out{1}.savedef.ofname = 'ea_inv_normparams.nii';
matlabbatch{1}.spm.util.defs.out{1}.savedef.savedir.saveusr = {directory};
spm_jobman('run',{matlabbatch});
