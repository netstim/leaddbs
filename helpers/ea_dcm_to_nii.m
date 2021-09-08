function ea_dcm_to_nii(method, dicom_dir, tmp_dir)

% simple wrapper that calls appropriate scripts to convert all dicoms in dicom_dir. Output is stored in tmp_dir
% method is an integer, 1 - 3
% 1 - dcm2nii, 2 - dicm2nii (Matlab), 3 - SPM)

switch method
    case 1 % dcm2niix
        ea_dcm2niix(dicom_dir, tmp_dir);
    case 2 % dicm2nii
        ea_dicm2nii(dicom_dir, tmp_dir);
    case 3 % SPM
        ea_spm_dicom_import(dicom_dir, tmp_dir);
end