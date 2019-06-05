function ea_delete(object, warn)
% wrapper for file/dir/graphic deleting, check existence beforehand

if nargin < 2
    warn = 0;
    warning('off');
end

if ~iscell(object)
    object = {object};
end

for i=1:numel(object)
    if exist(object{i}, 'file') == 2
        delete(object{i});
    elseif exist(object{i}, 'dir') == 7
        rmdir(object{i},'s');
    elseif contains(object{i}, '*') && ~isempty(dir(object{i}))
        contents = dir(object{i});
        for c=1:length(contents)
            fd = [contents(c).folder, filesep, contents(c).name];
            if exist(fd, 'file') == 2
                delete(fd);
            elseif exist(fd, 'dir') == 7
                rmdir(fd,'s');
            end
        end
    elseif warn
        warning([object{i}, ' does not exist!'])
    end
end
warning('on');
