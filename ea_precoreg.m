function ea_precoreg(fname,tname)
% pre-coreg predominant anat to MNI (without reslicing).
matlabbatch{1}.spm.spatial.coreg.estimate.ref = {[ea_space,tname,'.nii,1']};
matlabbatch{1}.spm.spatial.coreg.estimate.source = {[fname,',1']};
matlabbatch{1}.spm.spatial.coreg.estimate.other = {''};
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
spm_jobman('run',{matlabbatch});
clear matlabbatch