function parsedStruct = parseBIDSFilePath(filePath)
% Parse BIDS file path into a struct

if ~isBIDSFileName(filePath)
    error('Seems not a BIDS-like file name.')
else
    filePath = GetFullPath(filePath);
end

% Split file path into stripped path, file name and extension
strippedPath = regexp(filePath, ['.+(?=\', filesep, '[^\', filesep, ']+)'], 'match', 'once');
fullName = strrep(filePath, [strippedPath, filesep], '');
fileName = regexp(fullName, '[^.]+', 'match', 'once');
fileExt = strrep(fullName, fileName, '');

parsedStruct.dir = fileparts(strippedPath);

% Parse file name
entities = strsplit(fileName, '_');
if ~contains(entities{end}, '-')
    % filePath has the following patterns:
    % sub-XX_key1-value1_key2-value2_[modality].nii
    % sub-XX_key1-value1_key2-value2_[modality].nii.gz
    for i=1:length(entities)-1
        pair = regexp(entities{i}, '-', 'split', 'once');
        parsedStruct.(pair{1}) = pair{2};
    end
    parsedStruct.suffix = entities{end}; % Last one should be modality
else
    % filePath has the following patterns:
	% sub-XX_key1-value1_key2-value2.[ext]
    for i=1:length(entities)
        pair = regexp(entities{i}, '-', 'split', 'once');
        parsedStruct.(pair{1}) = pair{2};
    end
    parsedStruct.suffix = '';
end 

parsedStruct.ext = fileExt;
