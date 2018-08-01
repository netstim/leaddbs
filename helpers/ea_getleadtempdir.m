function tmpdir=ea_getleadtempdir

[~,uname]=fileparts(fileparts(ea_gethome));
tmpdir=[tempdir,uname,'_leaddbs',filesep];

if ~exist(tmpdir, 'dir')
    mkdir(tmpdir)
end
