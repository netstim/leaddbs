function ea_write_nii(nii)
nii.fname = GetFullPath(nii.fname);
spm_write_vol(nii,nii.img);
