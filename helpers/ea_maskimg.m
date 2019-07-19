function ea_maskimg(options,file,prefix)
directory=[options.root,options.patientname,filesep];
if ~exist([directory,'brainmask.nii'],'file')
	ea_genbrainmask(options);
end
[pth,fn,ext]=fileparts(file);
if ~exist([pth,filesep,prefix,fn,ext],'file')
    nii=ea_load_nii(file);
    bm=ea_load_nii([directory,'brainmask.nii']);
    nii.img=nii.img.*double(bm.img);
    nii.fname=[pth,filesep,prefix,fn,ext];
    spm_write_vol(nii,nii.img);
end
