function parsedStruct = parseBIDSFilePath(filePath)
% Parse BIDS file path into a struct

% Split file path into stripped path, file name and extension
[strippedPath, fileName, ext] = ea_niifileparts(GetFullPath(filePath));
parsedStruct.dir = fileparts(strippedPath);

% Parse file name
entities = strsplit(fileName, '_');
for i=1:length(entities)-1
    pair = regexp(entities{i}, '-', 'split', 'once');
    parsedStruct.(pair{1}) = pair{2};
end
parsedStruct.suffix = entities{end}; % Last one should be modality

parsedStruct.ext = ext;
