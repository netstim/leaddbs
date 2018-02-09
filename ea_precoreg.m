function ea_precoreg(fname,tname,options)
% pre-coreg predominant anat to MNI (without reslicing).
% set save_tmet = 1 if you want to store the initial MRI-transformation in
% the lead folder
save_tmat = 1;
%%
flags.cost_fun = 'nmi';
flags.sep = [4 2];
flags.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
flags.fwhm = [7 7];
tmp = spm_coreg([ea_space,tname,'.nii,1'],[fname,',1'],flags);
tmat = spm_matrix(tmp(:)');
if save_tmat
    save([options.root options.patientname filesep 'ea_precoregtransformation.mat'],'tmat');
end
anat_t1_vol = spm_vol(fname);
anat_t1_vol.mat = spm_get_space(fname, inv(tmat) * anat_t1_vol.mat);
spm_write_vol(anat_t1_vol,spm_read_vols(anat_t1_vol));