function ea_perform_lc(options)
% This function is the main execution function of LEAD-DBS Connectome.
% __________________________________________________________________________________
% Copyright (C) 2014 Charite University Medicine Berlin, Movement Disorders Unit
% Andreas Horn

%% structural parts
disp('*** Performing structural parts of LEAD-Connectome...');

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
    
    if ~exist([options.root,options.patientname,filesep,'connectomics'],'dir')
        mkdir([options.root,options.patientname,filesep,'connectomics']);
    end
    expfolder=[options.root,options.patientname,filesep,'connectomics',filesep,options.lc.general.parcellation,filesep];
    if ~exist(expfolder,'dir')
        mkdir(expfolder);
    end
    if ~exist([expfolder,'DTI_CM.mat'],'file')
        if ~exist([options.root,options.patientname,filesep,options.prefs.FTR_unnormalized],'file') % fibertracking has not been performed.
            warning('Fibertracking has not been done yet. Will do so before estimating structural connectivity matrix.');
            ea_perform_ft_proxy(options);
        end
        DTI_CM=ea_createCM_dti(options);
        cm=ea_export_CM_png(DTI_CM,'DTI Connectivity matrix',options,[0 10]);
        save([expfolder,'DTI_CM.mat'],'DTI_CM','-v7.3');
        saveas(cm,[expfolder,'DTI_CM.png']);
    end
end


disp('*** Done.');

%% functional parts

disp(['*** Performing functional parts of LEAD-Connectome...']);

% get files with rs-fMRI data
restfiles = dir([options.root,options.patientname,filesep,options.prefs.rest]);

% get number of files with rs-fMRI data
options.prefs.n_rest = numel(restfiles);

% display number of rs-fMRI files to analyze
disp(['*** ' num2str(options.prefs.n_rest) ' rs-fMRI files to analyze...']);

% the following steps clear existing connectome files which have been defined in `ea_prefs`
% because we will redefine the filenames in the code that follows
% ?? delete these names from `ea_prefs` eventually ??
options.prefs = rmfield(options.prefs,'rest');
options.prefs = rmfield(options.prefs,'pprest');
options.prefs = rmfield(options.prefs,'glrest');
options.prefs = rmfield(options.prefs,'gmtc');

% connectivity matrix steps:
if options.lc.func.compute_CM 
    if ~exist([options.root,options.patientname,filesep,'connectomics'],'dir')
        mkdir([options.root,options.patientname,filesep,'connectomics']);
    end
    % set export folder and check if it exists
    expfolder=[options.root,options.patientname,filesep,'connectomics',filesep,options.lc.general.parcellation,filesep];
    if ~exist(expfolder,'dir')
        mkdir(expfolder);
    end

    % set filenames for each rs-fMRI file
    for irest = 1:options.prefs.n_rest
        % set filenames for this iteration
        options.prefs.rest = restfiles(irest).name; 
        [~,name,ext] = fileparts(restfiles(irest).name); % get fileparts for other filenames
        options.prefs.pprest=strcat('sr',name,ext); % preprocessed rs-fMRI data
        options.prefs.glrest=strcat('glr',name,ext); % preprocessed and normalized rs-fMRI data
        options.prefs.gmtc=strcat(name,'_tc.mat'); % extracted timecourses of resting state fMRI data
        
        % create connectivity matrix for each rs-fMRI file
        if ~exist([expfolder,name,'_fMRI_CM.mat'],'file')
            disp(['Creating connectivity matrix for rs-fMRI file #',num2str(irest),': ',options.prefs.rest]);
            [fMRI_CM,gmtc]=ea_createCM_fmri(options);
            cm=ea_export_CM_png(fMRI_CM,['fMRI Connectivity matrix for ',name],options);
            save([expfolder,options.prefs.gmtc],'gmtc');
            save([expfolder,name,'_fMRI_CM.mat'],'fMRI_CM','-v7.3');
            saveas(cm,[expfolder,name,'_fMRI_CM.png']);
        end % end loop for this rs-fMRI file
    
    end % end loop for for all rs-fMRI files
    
end % end connectivity matrix section

disp(['*** Done analyzing ' num2str(options.prefs.n_rest) ' rs-fMRI files...']);

% graph metrics section below:
% for multiple file support: need to edit `ea_computeGM`
if options.lc.func.compute_GM || options.lc.struc.compute_GM % perform graph metrics
    if ~exist([options.root,options.patientname,filesep,'connectomics'],'dir')
        mkdir([options.root,options.patientname,filesep,'connectomics']);
    end
    expfolder=[options.root,options.patientname,filesep,'connectomics',filesep,options.lc.general.parcellation,filesep];
    if ~exist(expfolder,'dir')
        mkdir(expfolder);
    end
    ea_computeGM(options);
end

disp('*** Done.');






function ea_perform_ft_proxy(options)
%#function ea_ft_mesotracking_reisert
%#function ea_ft_deterministictracking_kroon
%#function ea_ft_globaltracking_reisert

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