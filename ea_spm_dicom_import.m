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

% Create temporary folder
tempFolder = fullfile(outdir, 'NIFTI_SPM');
if exist(tempFolder, 'dir')
    rmdir(tempFolder, 's');
end
mkdir(tempFolder);

if strcmp(mode, 'all')
    % Find all files under DICOM folder
    dcm = ea_regexpdir(dicomdir, '.*', 1, 'file');

    matlabbatch{1}.spm.util.import.dicom.data = dcm;
    matlabbatch{1}.spm.util.import.dicom.root = 'series';
    matlabbatch{1}.spm.util.import.dicom.outdir = {tempFolder};
    matlabbatch{1}.spm.util.import.dicom.protfilter = '.*';
    matlabbatch{1}.spm.util.import.dicom.convopts.format = 'nii';
    matlabbatch{1}.spm.util.import.dicom.convopts.meta = true;
    matlabbatch{1}.spm.util.import.dicom.convopts.icedims = 0;

    spm_jobman('run',{matlabbatch});
    clear matlabbatch
elseif strcmp(mode, 'serial')
    % Retrieve series of DCMs to be converted
    dcm = ea_searchdcm(dicomdir);

    for i=1:length(dcm)
        matlabbatch{1}.spm.util.import.dicom.data = dcm{i};
        matlabbatch{1}.spm.util.import.dicom.root = 'series';
        matlabbatch{1}.spm.util.import.dicom.outdir = {tempFolder};
        matlabbatch{1}.spm.util.import.dicom.protfilter = '.*';
        matlabbatch{1}.spm.util.import.dicom.convopts.format = 'nii';
        matlabbatch{1}.spm.util.import.dicom.convopts.meta = 0;
        matlabbatch{1}.spm.util.import.dicom.convopts.icedims = 0;

        spm_jobman('run',{matlabbatch});
        clear matlabbatch
    end
end

% Find all procotol folders
niiSubFolders = ea_regexpdir(tempFolder, '.*', 1, 'dir');

% Iterate through procotol folder, move NIfTI files to outdir
for i=1:numel(niiSubFolders)
    niiFiles = ea_regexpdir(niiSubFolders{i}, '\.nii$', 0);
    jsonFiles =  ea_regexpdir(niiSubFolders{i}, '\.json$', 0);
    [~, protocol] = fileparts(niiSubFolders{i});
    if numel(niiFiles) == 1
        movefile(niiFiles{1}, fullfile(outdir, [protocol, '.nii']));
        gzip(fullfile(outdir, [protocol, '.nii']));
        movefile(jsonFiles{1}, fullfile(outdir, [protocol, '.json']));
    elseif numel(niiFiles) > 1
        for f=1:numel(niiFiles)
            movefile(niiFiles{f}, fullfile(outdir, [protocol, '_', num2str(f, '%02d'), '.nii']));
            gzip(fullfile(outdir, [protocol, '_', num2str(f, '%02d'), '.nii']));
            movefile(jsonFiles{f}, fullfile(outdir, [protocol, '_', num2str(f, '%02d'), '.json']));
        end
    end
end

rmdir(tempFolder, 's');
