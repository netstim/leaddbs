function ea_spm_dicom_import(dicomdir, outdir)
% Wrapper to convert DICOM to NIfTI using SPM

if strcmp(dicomdir(end), filesep)
    dicomdir = dicomdir(1:end-1);
end

if nargin < 2
    outdir = fileparts(dicomdir);
end

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

% % Retrieve series of DCMs to be converted
% dcm = ea_searchdcm(dicomdir);

% for i=1:length(dcm)
%     matlabbatch{1}.spm.util.import.dicom.data = dcm{i};
%     matlabbatch{1}.spm.util.import.dicom.root = 'flat';
%     matlabbatch{1}.spm.util.import.dicom.outdir = {outdir};
%     matlabbatch{1}.spm.util.import.dicom.protfilter = '.*';
%     matlabbatch{1}.spm.util.import.dicom.convopts.format = 'nii';
%     matlabbatch{1}.spm.util.import.dicom.convopts.meta = 0;
%     matlabbatch{1}.spm.util.import.dicom.convopts.icedims = 0;
%
%     spm_jobman('run',{matlabbatch});
%     clear matlabbatch
% end
