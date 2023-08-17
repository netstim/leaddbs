function [fn,ext]=ea_niigz(base)
% wrapper for nifti file names. will return the filename present, if none
% present will return .nii (e.g. used for writing output nonexistent
% files)

[pth, name, ext]=ea_niifileparts(base);

nii = [pth, '.nii'];
niigz = [pth, '.nii.gz'];

if isfile(nii) && isfile(niigz)
    warning('off', 'backtrace');
    warning(['Duplicate .nii/.nii.gz files detected for ',name]);
    switch ext
        case '.nii' % explicitly asked for .nii
            fn = [pth,'.nii'];
            warning('Using .nii version since explicitly asked for.');
            ext = '.nii';
        otherwise
            fn = [pth,'.nii.gz'];
            ext = '.nii.gz';
    end
    warning('on', 'backtrace');
elseif ~isfile(nii) && isfile(niigz)
    fn = [pth,'.nii.gz'];
    ext = '.nii.gz';
elseif isfile(nii) && ~isfile(niigz)
    fn = [pth,'.nii'];
    ext = '.nii';
else % file not (yet) present, for now use .nii as default for output.
    fn = [pth,'.nii'];
    ext = '.nii';
end

