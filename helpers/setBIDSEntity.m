function filePath = setBIDSEntity(filePath, entity, value)
% Set BIDS entity value in file path

if strcmp(entity, 'suffix') % Set suffix
    filePath = regexprep(filePath, '_[a-zA-Z0-9]+(\.[a-zA-Z0-9\.]+$)', ['_', value, '$1']);
elseif strcmp(entity, 'ext') % Set extension
    filePath = regexprep(filePath, '_([a-zA-Z0-9]+\.)[a-zA-Z0-9\.]+$', ['_$1', value]);
elseif contains(filePath, [entity, '-']) % Entity already exist
    filePath = regexprep(filePath, [entity, '-[a-zA-Z0-9]+'], [entity, '-', value]);
else % Append new entity
    filePath = regexprep(filePath, '(_[a-zA-Z0-9]+\.[a-zA-Z0-9\.]+$)', ['_', entity, '-', value, '$1']);
end
