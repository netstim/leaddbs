function flag = isBIDSFileName(filePath)
% Check if file name is BIDS-like

[~, fileName] = fileparts(filePath);

if isempty(regexp(fileName, '^sub-.*_[a-zA-Z0-9]+(\.nii)?$', 'once'))
    flag = 0;
else
    flag = 1;
end
