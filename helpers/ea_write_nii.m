function ea_write_nii(nii)
nii.fname = GetFullPath(nii.fname);

% Fix endian in case missing
if numel(nii.dt) == 1
    [~, ~, endian] = computer;
    switch endian
        case 'L'
            nii.dt = [nii.dt, 0];
        case 'B'
            nii.dt = [nii.dt, 1];
    end
end

spm_write_vol(nii,nii.img);
