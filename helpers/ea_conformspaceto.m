function ea_conformspaceto(spacefn,toreslicefn,interp,mask)

sphdr=ea_open_vol(spacefn);
tohdr=ea_open_vol(toreslicefn);
if ~isequal(sphdr.mat,tohdr.mat) % volumes have different dimensions & hdr matrices.
    
    flags.mean=0;
    flags.which=1;
    flags.prefix='';
    if nargin>2
        flags.interp=interp;
    end
    if nargin>3
        flags.mask=mask;
    else
        flags.mask=0;
    end
    
    spm_reslice({spacefn,toreslicefn},flags);
    nii=ea_load_nii(toreslicefn);
    delete(toreslicefn);
    ea_write_nii(nii);
end
% make sure headers of images are exactly identical (also corrects for qform/sform issues).
sp=ea_load_untouch_nii(spacefn);
tr=ea_load_untouch_nii(toreslicefn);
sp.img=tr.img;
ea_save_untouch_nii(sp,toreslicefn);