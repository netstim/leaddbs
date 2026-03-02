function ea_write_nii(nii)
nii.fname = GetFullPath(nii.fname);
if endsWith(nii.fname, '.nii.gz')
    gzOutput = 1;
else
    gzOutput = 0;
end

% Ensure input to spm_write_vol has the ext of .nii
nii.fname = [ea_niifileparts(nii.fname), '.nii'];

% Fix endian in case missing
if isscalar(nii.dt)
    [~, ~, endian] = computer;
    switch endian
        case 'L'
            nii.dt = [nii.dt, 0];
        case 'B'
            nii.dt = [nii.dt, 1];
    end
end

spm_write_vol(nii,nii.img);

if gzOutput
    gzip(nii.fname);
    delete(nii.fname);
end

