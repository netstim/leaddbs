function ea_conformspaceto(spacefn,toreslicefn)

flags.mean=0;
flags.which=1;
flags.prefix='';

spm_reslice({spacefn,toreslicefn});