function ea_checkSpecialChars(paths)
% Check for special characters (and space) in the path

if ~isempty(paths) && strcmp(paths{1}, 'No Patient Selected')
    return
end

if isunix
    special_characters = cellfun(@(x) regexprep(x,'[\w-\/]',''), paths, 'uni', 0);
else
    special_characters = cellfun(@(x) regexprep(x(3:end),'[\w-\\]',''), paths, 'uni', 0);
end

warntxt = '';
for i = find(~cellfun(@isempty, special_characters))'
    warntxt = [warntxt sprintf('The folder: %s\ncontains the unsopported charaters (or space): ''%s''.\n\n', paths{i}, special_characters{i})];
end

if ~isempty(warntxt)
    warndlg([warntxt 'This can make some routines like normalization fail. Please rename the folders and select them again.'], 'Unsupported characters');
end
