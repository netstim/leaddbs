function ea_dicom_import(options)
% This function converts DICOM files in your in directory and outputs them
% to the out directory as specified by lead. This function is pretty much
% specialized to the DICOM format as used @ Charite University Medicine and
% might not work as expected in your clinical setting. You can
% alternatively import DICOM files using software like SPM or DCM2NII.
% __________________________________________________________________________________
% Copyright (C) 2014 Charite University Medicine Berlin, Movement Disorders Unit
% Andreas Horn

disp('Importing DICOM files...');

outdir = [options.root, options.patientname, filesep];

% check DICOM folder/zipfile under subject folder, can be named as:
% 'DICOM', 'DICOMDAT', 'DICOM.zip' or 'DICOMDAT.zip' (case insensitive).
dcmnames = ea_regexpdir(outdir, '^dicom(DAT)?(/|\\|\.zip)$', 0);

if isempty(dcmnames)
    % not found, suppose the subject folder is actually DICOM folder
    movefile([outdir, '*'],[outdir, 'DICOM'])
    try    movefile([outdir, 'DICOM', filesep, 'ea_ui.mat'], outdir); end % this isn't created when selecting multiple folders.
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

ea_dcm2niix(indir, outdir);

% delete DICOM folder
if options.prefs.dicom.dicomfiles
    rmdir(indir,'s');
end

% remove uncropped and untilted versions
fclean = ea_regexpdir(outdir, '(_Crop_1|_Tilt_1)\.nii$', 0);
fclean = unique(regexprep(fclean, '(_Crop_1|_Tilt_1)', ''));
for f=1:length(fclean)
    ea_delete(fclean{f});
end

if options.prefs.dicom.assign
    % assign image type here
    di = dir([outdir,'*.nii']);
    for d=1:length(di)
        dcfilename=[outdir,di(d).name];
        ea_imageclassifier({dcfilename});
    end
end


%% add methods dump:

ea_methods(options,['DICOM images were converted to the NIfTI file format, cropped and reoriented to standard NIfTI orientation using dcm2niiX software (e.g. see https://www.nitrc.org/projects/dcm2nii/).']);



