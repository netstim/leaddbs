function fn=ea_niigz(base)
% wrapper for nifti file names. will return the filename present

% rm base: 
[pth,base]=fileparts(base);


di=dir(fullfile(pth,[base,'*']));
if length(di)>1
        ea_error(['Duplicate .nii/.nii.gz files detected. Aborting. ',fullfile(pth,base),'.']);
end
[~,~,ext]=fileparts(di.name);
switch ext
    case '.gz'
        fn=fullfile(pth,[base,'.nii.gz']);
    case '.nii'
        fn=fullfile(pth,[base,'.nii']);
    otherwise
        fn='';
end
