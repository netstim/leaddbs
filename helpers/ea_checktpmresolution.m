function needstobebuilt=ea_checktpmresolution
options.earoot=ea_getearoot;
options.prefs=ea_prefs('');
needstobebuilt=1;
if ~exist([ea_space,'TPM.nii'],'file')
    return
end
% check resolution of TPM
V=ea_open_vol([ea_space,'TPM.nii']);
vox=ea_detvoxsize(V(1).mat);
if vox(1)==options.prefs.normalize.spm.resolution
    needstobebuilt=0;
end