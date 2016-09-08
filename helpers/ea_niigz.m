function fn=ea_niigz(base)
% wrapper for nifti file names. will return the filename present, if none
% present will return .nii.gz (e.g. used for writing output nonexistent
% files)

% rm base: 
[pth,base]=fileparts(base);
try
if strcmp(base(end-3:end),'.nii') % still has .nii ? .nii.gz has been applied
    [~,base]=fileparts(base);
end
end

nii=dir(fullfile(pth,[base,'.nii']));
niigz=dir(fullfile(pth,[base,'.nii.gz']));

if ~isempty(nii) && ~isempty(niigz)
    warning(['Duplicate .nii/.nii.gz files detected. ',fullfile(pth,base),'.']);
elseif isempty(nii) && ~isempty(niigz)
    fn=fullfile(pth,[base,'.nii.gz']);
elseif ~isempty(nii) && isempty(niigz)
    fn=fullfile(pth,[base,'.nii']);
else % file not (yet) present, for now use .nii as default for output.
    fn=fullfile(pth,[base,'.nii']);
end

