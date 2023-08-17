function ea_delete(object, warn)
% wrapper for file/dir deleting, check existence beforehand

if nargin < 2
    warn = 0;
    warning('off');
end

if ~iscell(object)
    object = {object};
end

for i=1:numel(object)
    if isfile(object{i})
        delete(object{i});
    elseif isfolder(object{i})
        rmdir(object{i},'s');
    elseif contains(object{i}, '*') && ~isempty(dir(object{i}))
        contents = dir(object{i});
        contents = contents(~ismember({contents.name}, {'.', '..'}));
        for c=1:length(contents)
            fd = [contents(c).folder, filesep, contents(c).name];
            if isfile(fd)
                delete(fd);
            elseif isfolder(fd)
                rmdir(fd,'s');
            end
        end
    elseif warn
        warning([object{i}, ' does not exist!'])
    end
end
warning('on');
