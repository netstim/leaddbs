function ea_conformspaceto(spacefn,toreslicefn,interp,mask,newfn)

sphdr=ea_open_vol(spacefn);
tohdr=ea_open_vol(toreslicefn);
if ~isequal(sphdr.mat,tohdr.mat) % volumes have different dimensions & hdr matrices.
    
    flags.mean=0;
    flags.which=1;
    if exist('newfn','var')
        flags.prefix='r';
    else
        flags.prefix='';
    end
    if nargin>2
        flags.interp=interp;
    end
    if nargin>3
        flags.mask=mask;
    else
        flags.mask=0;
    end
    
    spm_reslice({spacefn,toreslicefn},flags);
    if exist('newfn','var')
    [pth,fn,ext]=fileparts(toreslicefn);
    movefile(fullfile(pth,['r',fn,ext]),newfn);
    toreslicefn=newfn;
    end
    nii=ea_load_nii(toreslicefn);
    delete(toreslicefn);
    ea_write_nii(nii);
end
% make sure headers of images are exactly identical (also corrects for qform/sform issues).
sp=ea_load_untouch_nii(spacefn);
tr=ea_load_untouch_nii(toreslicefn);
sp.img=eval([class(tr.img),'(tr.img);']); % make sure to save data in same class as used before
sp.hdr.dime.bitpix=tr.hdr.dime.bitpix;
sp.hdr.dime.datatype=tr.hdr.dime.datatype; % keep datatype of original image.
ea_save_untouch_nii(sp,toreslicefn);