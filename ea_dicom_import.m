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




indir=[options.root,options.patientname,filesep,'DICOM',filesep];
outdir=[options.root,options.patientname,filesep];
tmpoutdir=outdir;

f=dir(indir);
% zipfile support..
for scan=1:length(f)
   [p,fi,e]=fileparts(f(scan).name);
   if strcmp(e,'.zip')
       unzip([indir,f(scan).name],indir);
       delete([indir,f(scan).name]);
   end
end


ea_dcm2niix(indir,outdir);

if options.prefs.dicom.dicomfiles % delete DICOM folder
    rmdir([options.root,options.patientname,filesep,'DICOM'],'s');
end


% remove uncropped versions
di=dir([outdir,'*_Crop_1.nii']);
for d=1:length(di)
    delete([outdir,di(d).name(1:end-11),'.nii']);    
end
% remove untilted versions
di=dir([outdir,'*_Tilt_1.nii']);
for d=1:length(di)
    delete([outdir,di(d).name(1:end-11),'.nii']);    
end

di=dir([outdir,'*.nii']);
for d=1:length(di)
    dcfilename=[outdir,di(d).name];
    ea_imageclassifier(dcfilename);
end


