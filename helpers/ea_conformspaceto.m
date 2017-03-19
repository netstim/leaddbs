function ea_conformspaceto(spacefn,toreslicefn,interp,mask)

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