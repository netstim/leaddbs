function parsedStruct = parseBIDSFilePath(filePath)
% Parse BIDS file path into a struct

if ~isBIDSFileName(filePath)
    error('Seems not a BIDS-like file name.')
else
    filePath = GetFullPath(filePath);
end

% Split file path into stripped path, file name and extension
% File name has the pattern of [\w-]+
% File extension has the pattern of (\.[^\W_]+){1,}$
parsedStruct.dir = regexp(filePath, ['.*(?=\', filesep, '[\w-]+(\.[^\W_]+){1,}$)'], 'match', 'once');
fileName = regexp(filePath, ['(?<=\', filesep, ')[\w-]+(?=(\.[^\W_]+){1,}$)'], 'match', 'once');
fileExt = regexp(filePath, '(\.[^\W_]+){1,}$', 'match', 'once');

% Parse file name
entities = strsplit(fileName, '_');
if ~contains(entities{end}, '-')
    % filePath has the following patterns (i.e., with suffix):
    % sub-XX_key1-value1_key2-value2_[modality].nii
    % sub-XX_key1-value1_key2-value2_[modality].nii.gz
    for i=1:length(entities)-1
        pair = regexp(entities{i}, '-', 'split', 'once');
        parsedStruct.(pair{1}) = pair{2};
    end
    parsedStruct.suffix = entities{end}; % Last one should be modality
else
    % filePath has the following patterns (i.e., without suffix):
	% sub-XX_key1-value1_key2-value2.[ext]
    for i=1:length(entities)
        pair = regexp(entities{i}, '-', 'split', 'once');
        parsedStruct.(pair{1}) = pair{2};
    end
    parsedStruct.suffix = '';
end 

parsedStruct.ext = fileExt;
