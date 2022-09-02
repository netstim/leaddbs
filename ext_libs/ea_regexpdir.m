function [list, isdir] = ea_regexpdir(rootdir, expstr, recursive, type)
% Wrapper for regexpdir (need to clear the persistent variable)

if ~exist('recursive','var')
    recursive = true;
end

% Search file or folder
if ~exist('type','var')
    type = 'file';
end

% Fix to use '^STR' pattern recursively
if recursive && startsWith(expstr, '^')
    if endsWith(expstr, '$')
        expstr = ['(^|.+[/\\])', expstr(2:end)];
    else
        expstr = ['(^|.+[/\\])', expstr(2:end), '[^/\\]*$'];
    end
end

clear regexpdir
list = regexpdir(rootdir, expstr, recursive);

switch type
    case {'f', 'file'}
        list = list(isfile(list));
        isdir = 0;
    case {'d', 'dir', 'folder'}
        % Remove filesep from the end of the folder path
        list = cellfun(@fileparts,list(isfolder(list)),'UniformOutput',false);
        isdir = 1;
    case {'a', 'all'}
        % Remove filesep from the end of the folder path
        list(isfolder(list)) = fileparts(list(isfolder(list)));
        isdir = zeros(size(list));
        isdir(isfolder(list)) = 1;
end

if ischar(list)
    list = {list};
end
