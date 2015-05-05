function ea_perform_lc(options)
% This function is the main execution function of LEAD-DBS Connectome.
% __________________________________________________________________________________
% Copyright (C) 2014 Charite University Medicine Berlin, Movement Disorders Unit
% Andreas Horn



%% structural parts


% perform fibertracking
if options.lc.struc.ft.do
    ea_perform_ft_proxy(options);
end

% normalize fibers
if options.lc.struc.ft.normalize % normalize fibertracts ? for now these should be denoted in Freiburg format.
    ea_normalize_fibers(options);
end


% create structural CM
if options.lc.struc.compute_CM
    mkdir([options.root,options.patientname,filesep,'connectomics']);
expfolder=[options.root,options.patientname,filesep,'connectomics',filesep,options.lc.general.parcellation,filesep];
mkdir(expfolder);
    if ~exist([options.root,options.patientname,filesep,options.prefs.FTR_unnormalized],'file') % fibertracking has not been performed.
    warning('Fibertracking has not been done yet. Will do so before estimating structural connectivity matrix.');
        ea_perform_ft_proxy(options);
    end  
    DTI_CM=ea_createCM_dti(options);
    cm=ea_export_CM_png(DTI_CM,'DTI Connectivity matrix',options);
    save([expfolder,'DTI_CM.mat'],'DTI_CM','-v7.3');
    saveas(cm,[expfolder,'DTI_CM.png']);
end




%% functional parts

if options.lc.func.compute_CM % create functional connectivity matrix
    mkdir([options.root,options.patientname,filesep,'connectomics']);
expfolder=[options.root,options.patientname,filesep,'connectomics',filesep,options.lc.general.parcellation,filesep];
mkdir(expfolder);
    [fMRI_CM,gmtc]=ea_createCM_fmri(options);
    cm=ea_export_CM_png(fMRI_CM,'fMRI Connectivity matrix',options);
    save([expfolder,options.prefs.gmtc],'gmtc');
    save([expfolder,'fMRI_CM.mat'],'fMRI_CM','-v7.3');
    saveas(cm,[expfolder,'fMRI_CM.png']);
end


if options.lc.func.compute_GM || options.lc.struc.compute_GM % perform graph metrics
    mkdir([options.root,options.patientname,filesep,'connectomics']);
expfolder=[options.root,options.patientname,filesep,'connectomics',filesep,options.lc.general.parcellation,filesep];
mkdir(expfolder);
    ea_computeGM(options);
end







function ea_perform_ft_proxy(options)
eval([options.lc.struc.ft.method,'(options)']); % triggers the fibertracking function and passes the options struct to it.
try load([options.root,options.patientname,filesep,'ea_ftmethod_applied']); end
if exist('ft_method_applied','var')
    try
        ft_method_applied{end+1}=options.lc.struc.ft.method;
    catch
        clear ft_method_applied
        ft_method_applied{1}=options.lc.struc.ft.method;
    end
else
    ft_method_applied{1}=options.lc.struc.ft.method;
end
ft_method_applied=options.lc.struc.ft.method;
save([options.root,options.patientname,filesep,'ea_ftmethod_applied'],'ft_method_applied');