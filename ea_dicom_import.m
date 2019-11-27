function ea_dicom_import(options)
% This function converts DICOM files in your in directory and outputs them
% to the out directory as specified by lead. This function is pretty much
% specialized to the DICOM format as used @ Charite University Medicine and
% might not work as expected in your clinical setting. You can
% alternatively import DICOM files using software like SPM or DCM2NII.
% __________________________________________________________________________________
% Copyright (C) 2014 Charite University Medicine Berlin, Movement Disorders Unit
% Andreas Horn

outdir = [options.root, options.patientname, filesep];

if isfield(options.dicomimp,'method') && strcmp(options.dicomimp.method,'BIDS nifti import (experimental)')
disp('Importing BIDS folder')

if ~exist('derivatives','dir')
    mkdir derivatives
end
    keyboard
    preopmri = dir('anat/*.nii.gz');

    
else
disp('Importing DICOM files...');

% check DICOM folder/zipfile under subject folder, can be named as:
% 'DICOM', 'DICOMDAT', 'DICOM.zip' or 'DICOMDAT.zip' (case insensitive).
dcmnames = ea_regexpdir(outdir, '^dicom(DAT)?(/|\\|\.zip)$', 0);

if isempty(dcmnames)
    % not found, suppose the subject folder is actually DICOM folder
    movefile([outdir, '*'],[outdir, 'DICOM'])
    try
        % this isn't created when selecting multiple folders.
        movefile([outdir, 'DICOM', filesep, 'ea_ui.mat'], outdir);
    end
    dcmname = [outdir, 'DICOM', filesep];
else % found DICOM folder/zipfile
    dcmname = dcmnames{1}; % only choose the first found one
end

if strcmp(dcmname(end-3:end), '.zip') % zip file under subject folder
    unzip(dcmname, [outdir, 'DICOM']);
    delete(dcmname);
    indir = [outdir, 'DICOM', filesep];
else % DICOM folder under subject folder
    zips = ea_regexpdir(dcmname, '.+\.zip$', 0); % check if DICOM folder contains zip files
    for i=1:numel(zips)
        unzip(zips{i}, dcmname);
        delete(zips{i});
    end
    indir = dcmname;
end

if isfield(options.dicomimp,'method')
    switch options.dicomimp.method
        case 1 % dcm2niix
            ea_dcm2niix(indir, outdir);
            ea_methods(options, ['DICOM images were converted to the '...
                'NIfTI file format using dcm2niix (see https://github.com/rordenlab/dcm2niix).']);
        case 2 % dicm2nii
            ea_dicm2nii(indir, outdir);
            ea_delete([outdir, 'dcmHeaders.mat']);
            ea_methods(options, ['DICOM images were converted to the '...
                'NIfTI file format using dicm2nii (see https://github.com/xiangruili/dicm2nii).']);
        case 3 % SPM
            ea_spm_dicom_import(indir, outdir);
            ea_methods(options, ['DICOM images were converted to the '...
                'NIfTI file format using SPM.']);
    end
else % use default set in prefs
    switch lower(options.prefs.dicom.tool)
        case 'dcm2niix'
            ea_dcm2niix(indir, outdir);
            ea_methods(options, ['DICOM images were converted to the '...
                'NIfTI file format using dcm2niix (see https://github.com/rordenlab/dcm2niix).']);
        case 'dicm2nii'
            ea_dicm2nii(indir, outdir);
            ea_delete([outdir, 'dcmHeaders.mat']);
            ea_methods(options, ['DICOM images were converted to the '...
                'NIfTI file format using dicm2nii (see https://github.com/xiangruili/dicm2nii).']);
        case 'spm' % SPM
            ea_spm_dicom_import(indir, outdir);
            ea_methods(options, ['DICOM images were converted to the '...
                'NIfTI file format using SPM.']);
    end
end

% delete DICOM folder
if options.prefs.dicom.dicomfiles
    rmdir(indir,'s');
end

% remove uncropped and untilted versions
fclean = ea_regexpdir(outdir, '(_Crop_1|_Tilt_1|_Tilt)\.nii$', 0);
fclean = unique(regexprep(fclean, '(_Crop_1|_Tilt_1)', ''));

for f=1:length(fclean)
    ea_delete(fclean{f});
end

end
