function fn=ea_niigz(base)
% wrapper for nifti file names. will return the filename present

% rm base: 
[pth,base]=fileparts(base);


if exist(fullfile(pth,[base,'.nii']),'file') && ~exist(fullfile(pth,[base,'.nii.gz']),'file')
    fn=fullfile(pth,[base,'.nii']);
elseif ~exist(fullfile(pth,[base,'.nii']),'file') && exist(fullfile(pth,[base,'.nii.gz']),'file')
    fn=fullfile(pth,[base,'.nii.gz']);
elseif exist(fullfile(pth,[base,'.nii']),'file') && exist(fullfile(pth,[base,'.nii.gz']),'file')
    ea_error(['Duplicate .nii/.nii.gz files detected. Aborting. ',fullfile(pth,base),'.']);
else % return empty string
    fn='';
end