function ea_coregmr_brainsfit(options)
% uses Brainsfit to coregister MRIs.

fprintf('\nInterpolating preoperative anatomical image\n')
ea_normalize_reslicepretra(options);
disp('Done.');

if exist([options.root,options.patientname,filesep,options.prefs.tranii_unnormalized],'file')
    disp('Coregistering postop MR tra to preop MRI...');
    ea_brainsfit([options.root,options.patientname,filesep,options.prefs.prenii_unnormalized],...
                 [options.root,options.patientname,filesep,options.prefs.tranii_unnormalized],...
                 [options.root,options.patientname,filesep,options.prefs.tranii_unnormalized],0);
    disp('Coregistration done.');
end

if exist([options.root,options.patientname,filesep,options.prefs.cornii_unnormalized],'file')
    disp('Coregistering postop MR tra to preop MRI...');
    ea_brainsfit([options.root,options.patientname,filesep,options.prefs.prenii_unnormalized],...
                 [options.root,options.patientname,filesep,options.prefs.cornii_unnormalized],...
                 [options.root,options.patientname,filesep,options.prefs.cornii_unnormalized],0);
    disp('Coregistration done.');
end

if exist([options.root,options.patientname,filesep,options.prefs.sagnii_unnormalized],'file')
    disp('Coregistering postop MR tra to preop MRI...');
    ea_brainsfit([options.root,options.patientname,filesep,options.prefs.prenii_unnormalized],...
                 [options.root,options.patientname,filesep,options.prefs.sagnii_unnormalized],...
                 [options.root,options.patientname,filesep,options.prefs.sagnii_unnormalized],0);
    disp('Coregistration done.');
end
