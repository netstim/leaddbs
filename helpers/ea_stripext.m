function fname = ea_stripext(fname)
% Strip file extension
% Note: for NIfTI image, '.nii' or '.nii.gz' will be stripped

fname = regexprep(fname, '(\.nii|\.nii\.gz|\.[^\.\s]+)$', '');
