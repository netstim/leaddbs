function ea_coreg_fa(options,usebrainmask)
% Coregister FA image to anatomical anchor image.
% Will export FA image from DTI if it doesn't exist yet.
%
% USAGE:
%
%    ea_coreg_fa(options,usebrainmask)
%
% INPUTS:
%    options:           LeadDBS options
%    usebrainmask:      bool, masking the coregistered image or not

directory=[options.root,options.patientname,filesep];
% check for presence of FA map
if ~exist([directory,options.prefs.fa2anat],'file')
    if ~exist([directory,options.prefs.fa],'file')
        if ~exist([directory,options.prefs.dti],'file')
            fprintf('\n\nNo dMRI data has been found. Proceeding without FA\n\n');
        else
            ea_isolate_fa(options);
        end
    end
    if exist([directory,options.prefs.fa],'file') % check again since could have been built above
        if exist([directory,options.prefs.fa],'file') % recheck if has been built.
            ea_backuprestore([directory,options.prefs.fa]);
            ea_coregimages(options,[directory,options.prefs.fa],[directory,options.prefs.prenii_unnormalized],[directory,options.prefs.fa2anat]);
        end
    end
end
if exist([directory,options.prefs.fa2anat],'file') % recheck if now is present.
    if usebrainmask
        ea_maskimg(options,[directory,options.prefs.fa2anat],bprfx);
    end
end
