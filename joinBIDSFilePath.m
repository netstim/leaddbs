function filePath = joinBIDSFilePath(parsedStruct)
% Join parsed BIDS file path struct (from parseBIDSFilePath) into file path

% Get dir
dir = parsedStruct.dir;

% Get key-value pairs
keys = fieldnames(parsedStruct);
keys = keys(2:end-2);
values = struct2cell(parsedStruct);
values = values(2:end-2);

% Construct file name
fileName = strjoin(strcat(keys, '-', values), '_');

% Get suffix
if isempty(parsedStruct.suffix)
    suffix = '';
else
    suffix = ['_', parsedStruct.suffix];
end

% Get extension
ext = parsedStruct.ext;

% Contruct file path
filePath = fullfile(dir, [fileName, suffix, ext]);
