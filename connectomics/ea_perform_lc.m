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
    if ~exist([expfolder,'DTI_CM.mat'],'file') || ...
       (isfield(options, 'overwriteapproved') && options.overwriteapproved)
        if ~exist([options.root,options.patientname,filesep,options.prefs.FTR_unnormalized],'file') % fibertracking has not been performed.
            warning('Fibertracking has not been done yet. Will do so before estimating structural connectivity matrix.');
            ea_perform_ft_proxy(options);
        end
        DTI_CM=ea_createCM_dti(options);
        cm=ea_export_CM_png(DTI_CM,'DTI Connectivity matrix',options,[0 mean(DTI_CM(:))+2*std(DTI_CM(:))]);
        save([expfolder,'DTI_CM.mat'],'DTI_CM','-v7.3');
        saveas(cm,[expfolder,'DTI_CM.png']);
    end
end

disp('*** Done.');

%% functional parts
if options.lc.func.compute_GM || options.lc.func.compute_CM
    disp('*** Performing functional parts of LEAD-Connectome...');

    % get files with rs-fMRI data
    restfiles = dir([options.root,options.patientname,filesep,options.prefs.rest_searchstring]);

    % get number of files with rs-fMRI data
    options.prefs.n_rest = numel(restfiles);

    % display number of rs-fMRI files to analyze
    disp(['*** ' num2str(options.prefs.n_rest) ' rs-fMRI files to analyze...']);
end

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
    [~, presentfiles] = ea_assignpretra(options);
    % set filenames for each rs-fMRI file
    for irest = 1:options.prefs.n_rest
        % set filenames for this iteration
        options.prefs.rest = restfiles(irest).name;
        [~,name,ext] = fileparts(restfiles(irest).name); % get fileparts for other filenames
        options.prefs.pprest=strcat('sr',name,ext); % preprocessed rs-fMRI data
        options.prefs.glrest=strcat('glr',name,ext); % preprocessed and normalized rs-fMRI data
        options.prefs.gmtc=strcat(name,'_tc.mat'); % extracted timecourses of resting state fMRI data

        % create connectivity matrix for each rs-fMRI file
        if ~ea_reglocked(options,['r',ea_stripext(options.prefs.rest),'_',ea_stripext(presentfiles{1})]) || ...
                ~exist([expfolder,ea_stripext(options.prefs.rest),'_fMRI_CM.mat'],'file');
            disp(['Creating connectivity matrix for rs-fMRI file #',num2str(irest),': ',options.prefs.rest]);
            [fMRI_CM, gmtc]=ea_createCM_fmri(options);
            cm=ea_export_CM_png(fMRI_CM,['fMRI Connectivity matrix for ',name],options);
            save([expfolder,options.prefs.gmtc],'gmtc');
            save([expfolder,name,'_fMRI_CM.mat'],'fMRI_CM','-v7.3');
            saveas(cm,[expfolder,name,'_fMRI_CM.png']);
        end % end loop for this rs-fMRI file
    end % end loop for for all rs-fMRI files
end % end connectivity matrix section

if options.lc.func.compute_GM || options.lc.func.compute_CM
    disp(['*** Done analyzing ' num2str(options.prefs.n_rest) ' rs-fMRI files...']);
end

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

    modes=cell(0);
    finas=cell(0);
    threshs=cell(0);
    fs=[]; % functional/structural measure (for structure-function similarity)
    if options.lc.func.compute_GM
        for irest = 1:options.prefs.n_rest
            modes{end+1}='fMRI';
            [~,rbase]=fileparts(restfiles(irest).name);
            finas{end+1}=[rbase,'_fMRI'];
            threshs{end+1}=options.lc.graph.fthresh;
            fs(end+1)=1;
        end
    end

    if options.lc.struc.compute_GM
        modes{end+1}='DTI';
        finas{end+1}='DTI';
        threshs{end+1}=options.lc.graph.sthresh;
        fs(end+1)=2;
    end

    try
        ea_computeGM(options,modes,finas,threshs,fs);
    catch
        ea_error('Please export connectivity matrices first.');
    end
end

disp('*** Done.');
