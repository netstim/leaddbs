function flag = isBIDSFileName(filePath)
% Check if file name is BIDS-like
%
% regexp will match:
%     sub-XX_key1-value1_key2-value2.[ext]
%     sub-XX_key1-value1_key2-value2_[modality].nii
%     sub-XX_key1-value1_key2-value2_[modality].nii.gz

[~, fileName, fileExt] = fileparts(filePath);
fullFileName = [fileName, fileExt];

pattern = '^sub-[^\W_]+(_[^\W_]+-[^\W_]+){1,}(_[^\W_]+)?(\.[^\W_]+){1,}$';

if isempty(regexp(fullFileName, pattern, 'once'))
    flag = false;
else
    flag = true;
end
