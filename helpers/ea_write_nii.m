function ea_write_nii(nii)
nii.fname = GetFullPath(nii.fname);

% ensure to output .nii no matter what was supplied (if supplying .nii.gz
% this leads to an error).
[pth,fn,ext]=fileparts(nii.fname);
nii.fname=fullfile(pth,[ea_stripext(fn),'.nii']);

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

if strcmp(ext,'.gz') % support for writing out .gz files
    gzip(nii.fname);
    delete(nii.fname);
end

