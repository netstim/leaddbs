function filePath = rmBIDSEntity(filePath, entity)
% Remove BIDS entity from file path

if ischar(entity)
    entity = {entity};
end

for i=1:length(entity)
    if strcmp(entity{i}, 'suffix') % Remove suffix
        filePath = regexprep(filePath, '_[a-zA-Z0-9]+(\.[a-zA-Z0-9\.]+$)', '$1');
    else % Remove other entities
        filePath = regexprep(filePath, [entity{i}, '-[a-zA-Z0-9]+_'], '');
    end
end
