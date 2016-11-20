function needstobebuilt=ea_checktpmresolution
options.earoot=ea_getearoot;
options.prefs=ea_prefs('');
needstobebuilt=1;
if ~exist([options.earoot,'templates',filesep,'TPM_2009b.nii'],'file')
    return
end
% check resolution of TPM
V=ea_open_vol([options.earoot,'templates',filesep,'TPM_2009b.nii']);
vox=ea_detvoxsize(V(1).mat);
if vox(1)==options.prefs.normalize.spm.resolution
    needstobebuilt=0;
end