function dirlist = ea_regexpdir(rootdir, expstr, recursive)
% wrapper for regexpdir (need to clear the persistent variable)

if ~exist('recursive','var')
    recursive = true;
end

% Fix to use '^STR' pattern recursively
if recursive && strcmp(expstr(1),'^')
    expstr = ['(^|.+[/\\])', expstr(2:end)];
end

clear regexpdir
dirlist = regexpdir(rootdir, expstr, recursive);
