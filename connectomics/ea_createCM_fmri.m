function [fMRI_CM,gmtc] = ea_createCM_fmri(options)

expfolder=[options.root,options.patientname,filesep,'connectomics',filesep,options.lc.general.parcellation,filesep];
% check if rest_preprocessing has been performed:
[~,restfname]=fileparts(options.prefs.rest);
[~,anatfname]=fileparts(options.prefs.prenii_unnormalized);
directory=[options.root,options.patientname,filesep];
if ~exist([directory,'sr',restfname,'.nii'],'file') ...
    || ~exist([directory,'r',restfname,'_c1',anatfname,'.nii'],'file') % preproecessing needs to be performed
    disp('No preprocessed fMRI-images found, processing...');
    ea_preprocess_fmri(options);
    disp('Done preprocessing fMRI data.');
end

if ~exist([expfolder,options.prefs.gmtc],'file')
    disp('No timecourses found, processing...');
    gmtc=ea_extract_timecourses(options);
else
    load([expfolder,options.prefs.gmtc]);
end

fMRI_CM=corrcoef(gmtc);
