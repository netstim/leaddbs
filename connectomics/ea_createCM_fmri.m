function [fMRI_CM,gmtc] = ea_createCM_fmri(options)

expfolder=[options.root,options.patientname,filesep,'connectomics',filesep,options.lc.general.parcellation,filesep];
% check if rest_preprocessing has been performed:
    if ~exist([options.root,options.patientname,filesep,'sr',options.prefs.rest],'file') % preproecessing needs to be performed
      disp('No preprocessed fMRI-images found, processing...');
        ea_preprocess_fmri(options); 
        disp('Done preprocessing fMRI data.');
    end

if ~exist([expfolder,options.prefs.gmtc],'file');
    disp('No timecourses found, processing...');
    gmtc=ea_extract_timecourses(options);
else
    load([expfolder,options.prefs.gmtc]);
end

fMRI_CM=corrcoef(gmtc,'rows','pairwise');




