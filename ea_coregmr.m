function ea_coregmr(options,automan)
% wrapper for coreg routines

% in CT imaging, coregistration is done elsewhere.
% also ignore when there is no tra/cor/sag existing (normal conn study)
if options.modality==2 || ~isfield(options.prefs,'tranii_unnormalized')
    return
end

directory=[options.root,options.patientname,filesep];

whichnormmethod=ea_whichnormmethod([options.root,options.patientname,filesep]);
if ismember(whichnormmethod,ea_getantsnormfuns)
    doreslice=1;
else
    doreslice=0;
end

% restore raw files -> postop files from prior attempts. & make backups
% from original files in any case.
try
if exist([directory,'raw_',options.prefs.tranii_unnormalized],'file')
    copyfile([directory,'raw_',options.prefs.tranii_unnormalized],[directory,options.prefs.tranii_unnormalized]);
else
    copyfile([directory,options.prefs.tranii_unnormalized],[directory,'raw_',options.prefs.tranii_unnormalized]);
end
end
try
    if exist([directory,'raw_',options.prefs.cornii_unnormalized],'file')
        copyfile([directory,'raw_',options.prefs.cornii_unnormalized],[directory,options.prefs.cornii_unnormalized]);
    else
        copyfile([directory,options.prefs.cornii_unnormalized],[directory,'raw_',options.prefs.cornii_unnormalized]);
    end
end
try
    if exist([directory,'raw_',options.prefs.sagnii_unnormalized],'file')
        copyfile([directory,'raw_',options.prefs.sagnii_unnormalized],[directory,options.prefs.sagnii_unnormalized]);
    else
        copyfile([directory,options.prefs.sagnii_unnormalized],[directory,'raw_',options.prefs.sagnii_unnormalized]);
    end
end

switch options.coregmr.method
    case 1 % SPM
        ea_coregmr_spm(options,automan,doreslice);
        return
    case 2 % ANTs
        ea_coregmr_ants(options);
        return
    case 3 % Brainsfit
        ea_coregmr_brainsfit(options);
        return
    case 4 % Hybrid SPM & ANTs
        ea_coregmr_spm(options,automan,0); % dont use doreslice here to refrain for doing two interpolations.
        ea_coregmr_ants(options);
    case 5 % Hybrid SPM & Brainsfit
        ea_coregmr_spm(options,automan,0); % dont use doreslice here to refrain for doing two interpolations.
        ea_coregmr_brainsfit(options);
    case 6 % Do nothing
        return
end
