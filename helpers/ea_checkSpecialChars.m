function found = ea_checkSpecialChars(paths)
% Check for special characters (and space) in the path

found = 0;

if isempty(paths)
    return
elseif ischar(paths)
    paths = {paths};
end

if iscolumn(paths)
    paths = paths';
end

if strcmp(paths{1}, 'No Patient Selected')
    return
end

if isunix
    special_characters = cellfun(@(x) regexprep(x,'[\w-\/]',''), paths, 'uni', 0);
else
    special_characters = cellfun(@(x) regexprep(x(3:end),'[\w-\\]',''), paths, 'uni', 0);
end

warntxt = '';
for i = find(~cellfun(@isempty, special_characters))
    warntxt = [warntxt sprintf('The folder: %s\ncontains the unsopported characters (or space): ''%s''.\n\n', paths{i}, special_characters{i})];
end

if ~isempty(warntxt)
    found = 1;
    warndlg([warntxt 'This is not recommended. Please rename and try again.'], 'Unsupported characters');
end
