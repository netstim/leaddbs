function tmpdir = ea_getleadtempdir

[~, uname] = fileparts(fileparts(ea_gethome));
tmpdir = [tempdir, uname, '_leaddbs', filesep];

ea_mkdir(tmpdir)
