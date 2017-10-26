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
    elseif warn
        warning([object{i}, ' does not exist!'])
    end
end
warning('on');
