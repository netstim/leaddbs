function list = ea_regexpdir(rootdir, expstr, recursive, type)
% Wrapper for regexpdir (need to clear the persistent variable)

if ~exist('recursive','var')
    recursive = true;
end

% Search file or folder
if ~exist('type','var')
    type = 'file';
end

% Fix to use '^STR' pattern recursively
if recursive && strcmp(expstr(1),'^')
    expstr = ['(^|.+[/\\])', expstr(2:end)];
end

clear regexpdir
list = regexpdir(rootdir, expstr, recursive);

switch type
    case {'f', 'file'}
        list = list(isfile(list));
    case {'d', 'dir', 'folder'}
        list = fileparts(list(isfolder(list)));
end

if ischar(list)
    list = {list};
end
