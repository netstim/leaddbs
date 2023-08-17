function [fMRI_CM,gmtc] = ea_createCM_fmri(options)

expfolder=[options.root,options.patientname,filesep,'connectomics',filesep,options.lc.general.parcellation,filesep];
% check if rest_preprocessing has been performed:
[~,restfname]=fileparts(options.prefs.rest);
[~,anatfname]=fileparts(options.prefs.prenii_unnormalized);
if ~ea_reglocked(options,['r',restfname,'_',anatfname]) || ...
        ~exist([options.root,options.patientname,filesep,'r',restfname,'_',anatfname,'.nii'],'file') % preprocessing needs to be performed
    disp('No (approved) preprocessed fMRI-images found, processing...');
    ea_preprocess_fmri(options);
    disp('Done preprocessing fMRI data.');
end

if ~ea_reglocked(options,['r',restfname,'_',anatfname]) || ~exist([expfolder,options.prefs.rest,'_tc.mat'],'file')

    disp('No timecourses found, processing...');
    gmtc=ea_extract_timecourses(options);
else
    load([expfolder,options.prefs.gmtc]);
end

fMRI_CM=corrcoef(gmtc);
