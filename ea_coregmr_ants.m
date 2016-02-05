function ea_coregmr_ants(options)
% uses ANTs instead of SPM to coregister MRIs.

disp('Interpolating preoperative anatomical image');
ea_normalize_reslicepretra(options);
disp('Done.');
disp('Coregistering postop MR tra to preop MRI...');
if exist([options.root,options.patientname,filesep,'raw_',options.prefs.tranii_unnormalized],'file')
copyfile([options.root,options.patientname,filesep,'raw_',options.prefs.tranii_unnormalized],[options.root,options.patientname,filesep,options.prefs.tranii_unnormalized]);
end
copyfile([options.root,options.patientname,filesep,options.prefs.tranii_unnormalized],[options.root,options.patientname,filesep,'raw_',options.prefs.tranii_unnormalized]);
ea_ants([options.root,options.patientname,filesep,options.prefs.prenii_unnormalized],...
          [options.root,options.patientname,filesep,options.prefs.tranii_unnormalized],...
          [options.root,options.patientname,filesep,options.prefs.tranii_unnormalized])
disp('Coregistration done.');

if exist([options.root,options.patientname,filesep,options.prefs.cornii_unnormalized],'file');
    disp('Coregistering postop MR tra to preop MRI...');
    if exist([options.root,options.patientname,filesep,'raw_',options.prefs.cornii_unnormalized],'file')
        copyfile([options.root,options.patientname,filesep,'raw_',options.prefs.cornii_unnormalized],[options.root,options.patientname,filesep,options.prefs.cornii_unnormalized]);
    end
    copyfile([options.root,options.patientname,filesep,options.prefs.cornii_unnormalized],[options.root,options.patientname,filesep,'raw_',options.prefs.cornii_unnormalized]);
    ea_ants([options.root,options.patientname,filesep,options.prefs.prenii_unnormalized],...
        [options.root,options.patientname,filesep,options.prefs.cornii_unnormalized],...
        [options.root,options.patientname,filesep,options.prefs.cornii_unnormalized])
    disp('Coregistration done.');
end

if exist([options.root,options.patientname,filesep,options.prefs.sagnii_unnormalized],'file');
    disp('Coregistering postop MR tra to preop MRI...');
    if exist([options.root,options.patientname,filesep,'raw_',options.prefs.sagnii_unnormalized],'file')
        copyfile([options.root,options.patientname,filesep,'raw_',options.prefs.sagnii_unnormalized],[options.root,options.patientname,filesep,options.prefs.sagnii_unnormalized]);
    end
        copyfile([options.root,options.patientname,filesep,options.prefs.sagnii_unnormalized],[options.root,options.patientname,filesep,'raw_',options.prefs.sagnii_unnormalized]);
    ea_ants([options.root,options.patientname,filesep,options.prefs.prenii_unnormalized],...
        [options.root,options.patientname,filesep,options.prefs.sagnii_unnormalized],...
        [options.root,options.patientname,filesep,options.prefs.sagnii_unnormalized])
    disp('Coregistration done.');
end