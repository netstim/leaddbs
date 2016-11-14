function [trimpath, basename, ext, slice] = ea_niifileparts(niifile)
% Return the nii file path and name without '.nii' or '.nii.gz' ext.
% Useful for FSL cli.
% For example, input '/PATH/TO/image.nii.gz,1' will return 
% ['/PATH/TO/image', 'image', 'nii.gz', ',1']

if regexp(niifile, '\.nii$', 'once') % 'image.nii'
    trimpath = niifile(1:end-4);
    ext = '.nii';
    slice = '';
elseif ~isempty(regexp(niifile, '(?<=\.nii),\d+$', 'match')) % 'image.nii,1'
    trimpath = niifile(1:regexp(niifile, '\.nii,\d+$', 'once')-1);
    ext = '.nii';
    slice = regexp(niifile, '(?<=\.nii),\d+$', 'match');
    slice = slice{:};
elseif regexp(niifile, '\.nii.gz$', 'once') % 'image.nii.gz'
    trimpath = niifile(1:end-7);
    ext = '.nii.gz';
    slice = '';
elseif ~isempty(regexp(niifile, '(?<=\.nii.gz),\d+$', 'match')) % 'image.nii.gz,1'
    trimpath = niifile(1:regexp(niifile, '\.nii.gz,\d+$', 'once')-1);
    ext = '.nii.gz';
    slice = regexp(niifile, '(?<=\.nii.gz),\d+$', 'match');
    slice = slice{:};
elseif ~isempty(regexp(niifile, ',\d+$', 'match')) % rare case: 'image,1'
    trimpath = niifile(1:regexp(niifile, ',\d+$', 'once')-1);
    ext = '';
    slice = regexp(niifile, ',\d+$', 'match');
    slice = slice{:};
else
    trimpath = niifile;
    ext = '';
    slice = '';
end

if isempty(fileparts(niifile))
    trimpath = ['.', filesep, trimpath];
end

[~, basename] = fileparts(trimpath);
