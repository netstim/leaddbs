function filePath = rmBIDSEntity(filePath, entity)
% Remove BIDS entity from file path

if ~isBIDSFileName(filePath)
    error('Seems not a BIDS-like file name.')
end

if ischar(entity)
    entity = {entity};
end

for i=1:length(entity)
    if strcmp(entity{i}, 'suffix') % Remove suffix
        filePath = regexprep(filePath, '_[^\W_]+((\.[^\W_]+){1,})$', '$1');
    else % Remove other entities
        filePath = regexprep(filePath, [entity{i}, '-[^\W_]+_'], '');
    end
end
