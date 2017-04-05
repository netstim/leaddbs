function ea_coregmr_ants(options,refine)
% uses ANTs to coregister MRIs.
if ~exist('refine','var')
    refine=0;
end

if exist([options.root,options.patientname,filesep,options.prefs.tranii_unnormalized],'file')
    disp('Coregistering postop MR tra to preop MRI...');
    ea_ants([options.root,options.patientname,filesep,options.prefs.prenii_unnormalized],...
            [options.root,options.patientname,filesep,options.prefs.tranii_unnormalized],...
            [options.root,options.patientname,filesep,options.prefs.tranii_unnormalized],...
            options.prefs.mrcoreg.writeoutcoreg,{},refine,options);
    disp('Coregistration done.');
end

if exist([options.root,options.patientname,filesep,options.prefs.cornii_unnormalized],'file')
    disp('Coregistering postop MR cor to preop MRI...');
    ea_ants([options.root,options.patientname,filesep,options.prefs.prenii_unnormalized],...
            [options.root,options.patientname,filesep,options.prefs.cornii_unnormalized],...
            [options.root,options.patientname,filesep,options.prefs.cornii_unnormalized],...
            options.prefs.mrcoreg.writeoutcoreg,{},refine,options);
    disp('Coregistration done.');
end

if exist([options.root,options.patientname,filesep,options.prefs.sagnii_unnormalized],'file')
    disp('Coregistering postop MR sag to preop MRI...');
    ea_ants([options.root,options.patientname,filesep,options.prefs.prenii_unnormalized],...
        [options.root,options.patientname,filesep,options.prefs.sagnii_unnormalized],...
        [options.root,options.patientname,filesep,options.prefs.sagnii_unnormalized],...
        options.prefs.mrcoreg.writeoutcoreg,{},refine,options);
    disp('Coregistration done.');
end
