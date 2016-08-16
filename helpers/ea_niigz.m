function fn=ea_niigz(base)
% wrapper for nifti file names. will return the filename present

% rm base: 
[pth,base]=fileparts(base);


nii=dir(fullfile(pth,[base,'.nii']));
niigz=dir(fullfile(pth,[base,'.nii.gz']));

if ~isempty(nii) && ~isempty(niigz)
    ea_error(['Duplicate .nii/.nii.gz files detected. Aborting. ',fullfile(pth,base),'.']);
elseif isempty(nii) && ~isempty(niigz)
    fn=fullfile(pth,[base,'.nii.gz']);
elseif ~isempty(nii) && isempty(niigz)
    fn=fullfile(pth,[base,'.nii']);
else % file not (yet) present, using .nii.gz as default
    fn=fullfile(pth,[base,'.nii.gz']);
end

