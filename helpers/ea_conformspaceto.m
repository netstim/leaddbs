function ea_conformspaceto(spacefn,toreslicefn,interp)

flags.mean=0;
flags.which=1;
flags.prefix='';
if nargin>2
    flags.interp=interp;
end

spm_reslice({spacefn,toreslicefn},flags);