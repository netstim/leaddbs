function tmpdir=ea_getleadtempdir

[~,uname]=fileparts(fileparts(ea_gethome));
tmpdir=[tempdir,uname,'_leaddbs',filesep];
warning('off')
mkdir(tmpdir)
warning('on');