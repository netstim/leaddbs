function needstobebuilt=ea_checktpmresolution
options.earoot=ea_getearoot;
options.prefs=ea_prefs('');
needstobebuilt=1;
if ~exist([ea_space,'TPM.nii'],'file')
    return
end

load([ea_space,'spacedef']);
if isfield(spacedef,'tpm')
   if strcmp(spacedef.tpm,'custom_fixed')
       needstobebuilt=0;
       return
   end
end

% check resolution of TPM
V = ea_open_vol([ea_space,'TPM.nii']);
if V.voxsize(1) == options.prefs.normalize.spm.resolution
    needstobebuilt=0;
end
