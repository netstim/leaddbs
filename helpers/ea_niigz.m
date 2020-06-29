function fn=ea_niigz(base)
% wrapper for nifti file names. will return the filename present, if none
% present will return .nii.gz (e.g. used for writing output nonexistent
% files)

[pth,name,ext]=ea_niifileparts(base);

nii=dir([pth,'.nii']);
niigz=dir([pth,'.nii.gz']);

if ~isempty(nii) && ~isempty(niigz)
    warning('off', 'backtrace');
    warning(['Duplicate .nii/.nii.gz files detected for ',name]);
    switch ext
        case '.nii' % explicitly asked for .nii
            fn=[pth,'.nii'];
            warning('Using .nii version since explicitly asked for.');
        otherwise
            fn=[pth,'.nii.gz'];
    end
    warning('on', 'backtrace');
elseif isempty(nii) && ~isempty(niigz)
    fn=[pth,'.nii.gz'];
elseif ~isempty(nii) && isempty(niigz)
    fn=[pth,'.nii'];
else % file not (yet) present, for now use .nii as default for output.
    fn=[pth,'.nii'];
end

