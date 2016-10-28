function ea_conv_wftr(options)
directory=[options.root,options.patientname,filesep];
if exist([directory,options.prefs.options.prefs.FTR_normalized],'file') 
    if ~exist([directory,'connectomes',filesep,'dMRI'],'file')
        mkdir([directory,'connectomes',filesep,'dMRI']);
    end
    movefile([directory,options.prefs.options.prefs.FTR_normalized],[directory,'connectomes',filesep,'dMRI',filesep,options.prefs.options.prefs.FTR_normalized]);
end