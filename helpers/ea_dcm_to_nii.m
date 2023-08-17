function niiFiles = ea_dcm_to_nii(dicom_dir, output_dir, tool)
% Simple wrapper that calls appropriate function to convert DICOM to NIfTI
% tool can be name of the tool or 1 - 4
% (1 - dcm2nii, 2 - dicm2nii, 3 - SPM, 4 - GDCM via Slicer)

if ~exist('output_dir', 'var') || isempty(output_dir)
    output_dir = fullfile(dicom_dir, 'tmp');
end

ea_mkdir(output_dir);

% catch the case where method is a string
if isnumeric(tool)
    switch tool
        case 1
            tool = 'dcm2niix';
        case 2
            tool = 'dicm2nii';
        case 3
            tool = 'spm';
        case 4
            tool = 'gdcm via slicer';
        otherwise
            tool = 'dcm2niix';
    end
else
    if ~ismember(lower(tool), {'dcm2niix', 'dicm2nii', 'spm', 'gdcm via slicer'})
        tool = 'dcm2niix';
    end
end

switch lower(tool)
    case 'dcm2niix'
        ea_dcm2niix(dicom_dir, output_dir);
    case 'dicm2nii'
        ea_dicm2nii(dicom_dir, output_dir);
    case 'spm'
        ea_spm_dicom_import(dicom_dir, output_dir);
    case 'gdcm via slicer'
        ea_gdcm_slicer_dicom_import(dicom_dir, output_dir);
end

niiFiles = ea_regexpdir(output_dir, '\.nii\.gz$', 0);
