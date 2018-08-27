function ea_spm_dicom_import(dicomdir, outdir, mode)
% Wrapper to convert DICOM to NIfTI using SPM
%
% 'mode' can be either 'all' (DEFAULT) or 'serial'. If it's set to 'all',
% all the DCMs found within the 'dicomdir' will be feed to the SPM batch
% job. If it's set to 'serial', DCMs will be sorted for different series
% first (using function from dicm2nii by Xiangrui Li), the conversion will
% then be performed series by series.

if strcmp(dicomdir(end), filesep)
    dicomdir = dicomdir(1:end-1);
end

if nargin < 2
    outdir = fileparts(dicomdir);
end

if nargin < 3
    mode = 'all'; % Feed all DCMs to SPM batch job
end

if strcmp(mode, 'all')
    % Find all files under DICOM folder
    dcm = ea_regexpdir(dicomdir, ['[^\', filesep, ']$'], 1);

    matlabbatch{1}.spm.util.import.dicom.data = dcm;
    matlabbatch{1}.spm.util.import.dicom.root = 'flat';
    matlabbatch{1}.spm.util.import.dicom.outdir = {outdir};
    matlabbatch{1}.spm.util.import.dicom.protfilter = '.*';
    matlabbatch{1}.spm.util.import.dicom.convopts.format = 'nii';
    matlabbatch{1}.spm.util.import.dicom.convopts.meta = 0;
    matlabbatch{1}.spm.util.import.dicom.convopts.icedims = 0;

    spm_jobman('run',{matlabbatch});
    clear matlabbatch
elseif strcmp(mode, 'serial')
    % Retrieve series of DCMs to be converted
    dcm = ea_searchdcm(dicomdir);

    for i=1:length(dcm)
        matlabbatch{1}.spm.util.import.dicom.data = dcm{i};
        matlabbatch{1}.spm.util.import.dicom.root = 'flat';
        matlabbatch{1}.spm.util.import.dicom.outdir = {outdir};
        matlabbatch{1}.spm.util.import.dicom.protfilter = '.*';
        matlabbatch{1}.spm.util.import.dicom.convopts.format = 'nii';
        matlabbatch{1}.spm.util.import.dicom.convopts.meta = 0;
        matlabbatch{1}.spm.util.import.dicom.convopts.icedims = 0;

        spm_jobman('run',{matlabbatch});
        clear matlabbatch
    end
end
