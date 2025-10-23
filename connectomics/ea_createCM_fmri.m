function [fMRI_CM,gmtc] = ea_createCM_fmri(options)

expfolder=[options.root,options.patientname,filesep,'connectomics',filesep,options.lc.general.parcellation,filesep];
% check if rest_preprocessing has been performed:
[~,restfname]=fileparts(options.prefs.rest);
[~,anatfname]=fileparts(options.prefs.prenii_unnormalized);

% BIDS FIX: Use simple file check for preprocessed fMRI
% For BIDS paths, the preprocessed files are in the same directory as the raw file
directory = [options.root, options.patientname, filesep];

% Build correct paths using ea_prependFilename
rest_sr = ea_prependFilename(options.prefs.rest, 'sr');
rest_r = ea_prependFilename(options.prefs.rest, 'r');

preprocFile1 = fullfile(directory, rest_sr);
preprocFile2 = fullfile(directory, rest_r);

if ~exist(preprocFile1, 'file') && ~exist(preprocFile2, 'file')
    disp('No (approved) preprocessed fMRI-images found, processing...');
    ea_preprocess_fmri(options);
    disp('Done preprocessing fMRI data.');
else
    disp('Preprocessed fMRI files found, skipping preprocessing.');
end

% Check if timecourses exist
tcFile = [expfolder, options.prefs.gmtc];
if ~exist(tcFile, 'file')
    disp('No timecourses found, processing...');
    gmtc=ea_extract_timecourses(options);
else
    disp('Loading existing timecourses...');
    load(tcFile);
end

fMRI_CM=corrcoef(gmtc);
