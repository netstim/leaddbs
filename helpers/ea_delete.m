function ea_delete(filename)
% wrapper for file deleting, check existence beforehand

if ~iscell(filename)
    filename = {filename};
end

for i=1:numel(filename)
    if exist(filename{i}, 'file')
        delete(filename{i});
    end
end
