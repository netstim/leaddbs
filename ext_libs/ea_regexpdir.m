function dirlist = ea_regexpdir(rootdir, expstr, recursive)
% wrapper for regexpdir (need to clear the persistent variable)

clear regexpdir
dirlist = regexpdir(rootdir, expstr, recursive);
