function fname = ea_stripext(fpath)
% Extract file name from path, strip file extension
% Note: for NIfTI image, '.nii' or '.nii.gz' will be stripped

fname = regexp(fpath, ['([^\', filesep, '\f\n\r\t\v]+?)(?=(\.nii(\.gz)?|\.[^\.\s]+)?$)'], 'match', 'once');
