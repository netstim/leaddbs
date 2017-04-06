function ea_resliceanat(options)
% Reslice the preoperative anatomical image if the resolution is not sufficient

directory=[options.root,options.patientname,filesep];
V=spm_vol([directory,options.prefs.prenii_unnormalized]);

dim=V.mat(logical(eye(4)));
dim=abs(dim(1:3));
if any(dim>0.7)
    fprintf('\nInterpolating preoperative anatomical image...')
    ea_reslice_nii([directory,options.prefs.prenii_unnormalized],[directory,options.prefs.prenii_unnormalized],[0.7 0.7 0.7],0,[],1);
    fprintf('\n\n');
end
