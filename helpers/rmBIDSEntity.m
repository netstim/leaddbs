function filePath = rmBIDSEntity(filePath, entity)
% Remove BIDS entity from file path

if ~isBIDSFileName(filePath)
    error('Seems not a BIDS-like file name.')
end

if ischar(entity)
    entity = {entity};
end

% Parse file path for further operation
parsedStruct = parseBIDSFilePath(filePath);

for i=1:length(entity)
    if strcmp(entity{i}, 'suffix') % Remove suffix
        parsedStruct.suffix = '';
    else % Remove other entities
        if isfield(parsedStruct, entity{i})
            parsedStruct = rmfield(parsedStruct, entity{i});
        else
            warning('off', 'backtrace');
            warning('Entity ''%s'' does not exist!', entity{i});
            warning('on', 'backtrace');
        end
    end
end

% Join BIDS file path struct into file path
filePath = joinBIDSFilePath(parsedStruct);
