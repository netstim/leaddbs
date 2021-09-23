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

% Parse file path for further operation
parsedStruct = parseBIDSFilePath(filePath);

for i=1:length(entities)
    entity = entities{i};
    value = values{i};
    if strcmp(entity, 'suffix') % Set suffix
        parsedStruct.suffix = value;
    elseif strcmp(entity, 'ext') % Set extension
        if ~startsWith(value, '.') % 'ext' instead of '.ext' provided
            parsedStruct.ext = ['.', value];
        else
            parsedStruct.ext = value;
        end
    elseif isfield(parsedStruct, entity) % Entity already exist
        parsedStruct.(entity) = value;
    else % Append new entity
        % Insert key-value pair before suffix
        parsedEntities = fieldnames(parsedStruct);
        parsedEntities = [parsedEntities(1:end-2); entity; parsedEntities(end-1:end)];
        parsedValues = struct2cell(parsedStruct);
        parsedValues = [parsedValues(1:end-2); value; parsedValues(end-1:end)];

        % Construct new parsed struct
        parsedStruct = cell2struct(parsedValues, parsedEntities, 1);
    end
end

% Join BIDS file path struct into file path
filePath = joinBIDSFilePath(parsedStruct);
