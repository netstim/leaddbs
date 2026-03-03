function ea_delete(object, opts)
% Enhanced file/folder deletion function

arguments
    object % Files/folders to be deleted
    opts.warn {mustBeNumericOrLogical} = 0 % Show warning or not
    opts.recycle {mustBeMember(opts.recycle, {'', 'off', 'on'})} = '' % Use recycle or not
end

% Update the recycle status
previousState = recycle;
if ~isempty(opts.recycle)
    recycle(opts.recycle);
end
    
if ~iscell(object)
    object = {object};
end

for i = 1:numel(object)
    obj = object{i};
    if ~(ischar(obj) || isstring(obj))
        continue;
    else
        obj = ea_strip(char(obj));
    end

    % Global safety-guard: do not delete empty or pure wildcards
    % without path definitions
    if isempty(obj) || strcmp(obj, '*')
        if opts.warn
            ea_cprintf('CmdWinWarnings', 'ea_delete: skip deleting "%s".\n', obj);
        end
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
    elseif opts.warn
        ea_cprintf('CmdWinWarnings', '"%s" does not exist!\n', obj);
    end
end

% Restore recycle status in user preference
if ~isempty(opts.recycle)
    recycle(previousState);
end
