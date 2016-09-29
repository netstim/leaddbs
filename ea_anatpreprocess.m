function ea_anatpreprocess(fpth)
[pth,fn,ext]=fileparts(fpth);
copyfile(fpth,fullfile(pth,['k',fn,ext]));
ea_dcm2nii(fpth);

V=ea_open_vol(fpth);
if any(V.dim<20)
    ea_warning('DCM2NII based autocrop and -reorient failed. Continuing without.');
    copyfile(fullfile(pth,['k',fn,ext]),fpth);
end
ea_bias_field_correction(fpth)
V=ea_open_vol(fpth);
if any(V.dim<20)
    ea_warning('Biasfield correction failed. Continuing without.');
    copyfile(fullfile(pth,['k',fn,ext]),fpth);
end
delete(['k',fn,ext]);

