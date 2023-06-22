function [list, isdir] = ea_regexpdir(rootdir, expstr, recursive, type, listHidden)
% Wrapper for regexpdir (need to clear the persistent variable)

% Recursive search by default
if ~exist('recursive', 'var')
    recursive = true;
end

% Search file or folder
if ~exist('type', 'var')
    type = 'file';
end

% Do not list hidden files by default
if ~exist('listHidden', 'var')
    listHidden = false;
end

% Rewrite pattern when '^STR$' used
if startsWith(expstr, '^')
    if endsWith(expstr, '$')
        expstr = ['(^|.+\', filesep, ')', expstr(2:end-1), '\', filesep, '?$'];
    else
        expstr = ['(^|.+\', filesep, ')', expstr(2:end)];
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
        list = fileparts(list(isfolder(list)));
        isdir = 1;
    case {'a', 'all'}
        % Remove filesep from the end of the folder path
        list(isfolder(list)) = erase(list(isfolder(list)), filesep + textBoundary('end'));
        isdir = zeros(size(list));
        isdir(isfolder(list)) = 1;
end

if ischar(list)
    list = {list};
end

if ~listHidden
    [~, fname] = fileparts(list);
    list(startsWith(fname, '.')) = [];
    if ismember(type, {'a', 'all'})
        isdir(startsWith(fname, '.')) = [];
        if isempty(isdir)
            isdir = 0;
        end
    end
end
