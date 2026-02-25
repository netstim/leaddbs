function ea_delete(object, warn)
% wrapper for file/dir deleting, check existence beforehand

if nargin < 2
    warn = 0;
    warning('off');
end

if ~iscell(object)
    object = {object};
end

for i = 1:numel(object)
    obj = object{i};
    if ~(ischar(obj) || isstring(obj))
        continue;
    else
        obj = strip(char(obj));
    end

    % Global safety-guard: do not delete empty or pure wildcards
    % without path definitions
    if isempty(obj) || strcmp(obj, '*')
        % optional: warning('ea_delete: refusing to delete "%s".', obj);
        continue;
    end

    % go on
    if isfile(obj)
        delete(obj);
    elseif isfolder(obj)
        rmdir(obj, 's');
    elseif contains(obj, '*') && ~isempty(dir(obj))
        contents = dir(obj);
        contents = contents(~ismember({contents.name}, {'.', '..'}));
        for c = 1:length(contents)
            fd = [contents(c).folder, filesep, contents(c).name];
            if isfile(fd)
                delete(fd);
            elseif isfolder(fd)
                rmdir(fd, 's');
            end
        end
    elseif warn
        warning([obj, ' does not exist!'])
    end
end

warning('on');
