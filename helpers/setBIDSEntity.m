function filePath = setBIDSEntity(filePath, varargin)
% Set BIDS entity and value in file path

if ~isBIDSFileName(filePath)
    error('Seems not a BIDS-like file name.')
end

entities = varargin(1:2:end-1);
values = varargin(2:2:end);

if length(entities) ~= length(values)
    error('Length of entities doesn''t match length of values!')
end

for i=1:length(entities)
    entity = entities{i};
    value = values{i};
    if strcmp(entity, 'suffix') % Set suffix
        filePath = regexprep(filePath, '_[^\W_]+(\.[a-zA-Z0-9\.]+$)', ['_', value, '$1']);
    elseif strcmp(entity, 'ext') % Set extension
        filePath = regexprep(filePath, '_([^\W_]+\.)[a-zA-Z0-9\.]+$', ['_$1', value]);
    elseif contains(filePath, [entity, '-']) % Entity already exist
        filePath = regexprep(filePath, [entity, '-[^\W_]+'], [entity, '-', value]);
    else % Append new entity
        filePath = regexprep(filePath, '(_[^\W_]+\.[a-zA-Z0-9\.]+$)', ['_', entity, '-', value, '$1']);
    end
end
