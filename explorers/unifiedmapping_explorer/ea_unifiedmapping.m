classdef ea_unifiedmapping < handle
    % Discriminative fiber class to handle visualizations of discriminative fibers in lead dbs resultfig / 3D Matlab figures
    % A. Horn

    properties (SetObservable)
        tool = 1; %option of 1 = sweetspotmapping, 2 = fiberfiltering, 3 = networkmapping
        fileformatversion; % 1.2 is current format
        M % content of lead group project
        resultfig % figure handle to plot results
        ID % name / ID of unified analysis object
        posvisible = 1 % pos tract/vox visible
        negvisible = 0 % neg tract/vox visible
        roivisible = 0 % show ROI (usually VTAs)
        connfibvisible = 0 % show all connected tracts in white
        showposamount = [25 25] % two entries for right and left
        shownegamount = [25 25] % two entries for right and left
        statmetric = nan;
        statsettings = struct;
        calcsettings = struct; %calcthreshold, connectivity_type
        threshstrategy = 'Percentage Relative to Peak'; % can be 'Relative to Amount' or 'Fixed Amount'
        multi_pathways = 0 % if structural connectome is devided into pathways (multiple .mat in dMRI_MultiTract)
        map_list % list that contains global indices of the first fibers in each pathway (relevant when multi_pathways = 1)
        pathway_list % list that contains names of pathways (relevant when multi_pathways = 1)
        connFiberInd % legacy
        switch_connectivity = 0 % flag if connectivity type was changed in the GUI
        nestedLOO = false       % if true, will conducted LOO in the training set
        corrtype = 'Spearman' % correlation strategy in case of statmetric ==  Correlations / E-fields (Irmen 2020). In case of one-sample tests used for 'T-Tests' vs 'Wicoxon Signed Rank Tests'.
        SigmoidTransform = 0;
        multitractmode = 'Single Tract Analysis' % multi mode now coded by this value %should we use abreviation?
        numpcs = 4; % standard value of how many PCs to compute in case of PCA mode
        doactualprediction = 0; % set up nested CVs to carry out actual predictions of response variables
        predictionmodel = 'Linear'; % type of glm used to fit fiber values to actual scores
        showsignificantonly = 0
        alphalevel = 0.05
        multcompstrategy = 'FDR'; % could be 'Bonferroni'
        subscore
        explorerdrawn
        results = struct
        customRoi = struct % struct used only for pseudoM case (customRoi.isbinary and customRoi.minmax)
        activateby={}; % entry to use to show fiber activations
        cvlivevisualize = 0; % if set to 1 shows crossvalidation results during processing.
        basepredictionon = 'Mean of Scores';
        fiberdrawn % struct contains fibercell and vals drawn in the resultfig
        spotdrawn % struct contains nifti of the sweetspot
        networkdrawn %struct contains nifti of the networkmap
        vizmode
        smooth_fp = 0; %for networkmapping, smooth fingerprints
        normalize_fp = 0; %for networkmapping, normalize fingerprints
        mask_vta_fp = 0; %for networkmapping, remove VTA voxels from fingerprints
        cvmask = 'Gray Matter';
        model='Smoothed'; %for networkmapping
        drawobject = struct % actual streamtube handle
        drawvals % weights of the fibers drawn
        connfiberdrawn % struct contains white connected fibers
        conndrawobject % actial streamtube handle for the latter
        roidrawobject % actual patch handle for ROI/VTAs
        roidata % data used for ROIs (nifti file cell usually w/ Efields
        roiprotocol % protocol for drawing rois used to check if need to be redrawn.
        patientselection % selected patients to include. Note that connected fibers are always sampled from all (& mirrored) VTAs of the lead group file
        setlabels={};
        setselections={};
        customselection % selected patients in the custom test list
        allpatients % list of all patients (as from M.patient.list)
        mirrorsides = 0 % flag to mirror VTAs / Efields to contralateral sides using ea_flip_lr_nonlinear()
        responsevar % response variable
        responsevarlabel % label of response variable
        covars = {} % covariates
        covarlabels = {} % covariate labels
        analysispath % where to store results
        leadgroup % redundancy protocol only, path to original lead group project
        useExternalModel = false
        ExternalModelFile = 'None'
        
        % NM visualization
        NMviz=struct; 
        % NMviz.vizmode='Regions'; % way to plot results
        % NMviz.model='Smoothed'; % in case of surface above, on which surface to plot.
        % NMviz.modelLH=1; % show left hemisphere
        % NMviz.modelRH=1; % show right hemisphere


        % stats: (how many fibers available and shown etc for GUI)
        modelNormalization = 'None';
        numBins=15;
        stats
        
        % additional settings:
        rngseed = 'default';
        Nperm = 1000 % how many permutations in leave-nothing-out permtest strategy
        kfold = 5 % divide into k sets when doing k-fold CV
        Nsets = 5 % divide into N sets when doing Custom (random) set test
        adjustforgroups = 1 % adjust correlations for group effects
        kIter = 1;
        roiintersectdata = {}; %roi, usually efield with which you can calculate fiber intersection
        roithresh = 200; % threshold above which efield metrics are considered
        drawTool = 'sweetspotmapping'; % active draw tool: sweetspotmapping, fiberfiltering, networkmapping
        activated = struct; % activated.fiberfiltering, activated.sweetspotmapping, activated.networkmapping % for activation status mapping.
        
        % misc
        runwhite = 0; % flag to calculate connected tracts instead of stat tracts
        e_field_metric = 'Magnitude'; % 'Magnitude' or 'Projection'

        %color
        colorbar % colorbar information
        posBaseColor = [1,1,1] % positive main color
        poscolor = [0.9176,0.2000,0.1373] % positive peak color
        negBaseColor = [1,1,1] % negative main color
        negcolor = [0.2824,0.6157,0.9725] % negative peak color
        hasResults = false; % results check
        AdditionalSettingsSavePath = [];
      
    end 

    properties (Access = private)
        switchedFromSpace=3 % if switching space, this will protocol where from
    end

    methods
        function obj=ea_unifiedmapping(analysispath) % class constructor
            if exist('analysispath', 'var') && ~isempty(analysispath)
                obj.analysispath = analysispath;
                [~, ID] = fileparts(obj.analysispath);
                obj.ID = ID;
            end
        end

        function initialize(obj,datapath,resultfig)
            % statsettings
            % initial hard threshold to impose on (absolute) nifti files only when calculating the data
            obj.statsettings.doVoxels = 1;
            obj.statsettings.doFibers = 1;
            obj.statsettings.outcometype = 'gradual';
            obj.statsettings.stimulationmodel = 'Electric Field';
            obj.statsettings.efieldmetric = 'Sum'; % if statmetric == ;Correlations / E-fields (Irmen 2020)’, efieldmetric can calculate sum, mean or peak along tracts
            obj.statsettings.efieldthreshold_spot = 200;
            obj.statsettings.efieldthreshold_tract = 200;
            obj.statsettings.efieldthreshold_network = 50;
            % Optional legacy sweetspot behavior. When empty, subthreshold
            % values remain explicit zeros rather than being treated as
            % missing observations.
            obj.statsettings.nanthreshold = [];
            obj.statsettings.sweetspotresolution = 0.5; % resolution of sweetspot in mm
            obj.statsettings.connthreshold = 20;
            obj.statsettings.statfamily = 'Correlations'; % the
            obj.statsettings.stattest = 'Spearman';
            obj.statsettings.H0 = 'Average';
            obj.calcsettings.selectedTool = 1; %1 = sweetspotmapping, 2 = fiberfiltering, 3 = networkmapping;
            obj.calcsettings.calcthreshold = 50;
            obj.calcsettings.switch_connectivity = 1;
            obj.calcsettings.connectivity_type = 1; %1 = vta, 2 = PAM
            obj.calcsettings.functionalresolution = '2 mm'; 
            obj.calcsettings.structuralresolution = '2 mm';
            obj.calcsettings.calcmethod = 1; %1 = e-field based method, 2 = fiber based method
            obj.calcsettings.calcspace = 1; %0 = native space, 1 = MNI space
            obj.calcsettings.netmap_connectome = '';
            obj.calcsettings.fibfilt_connectome = '';
            if ~isfield(obj.NMviz,'modelRH')
            obj.NMviz.modelRH=1;
            end
            if ~isfield(obj.NMviz,'modelLH')
            obj.NMviz.modelLH=1;
            end
            obj.subscore.vars = {};
            obj.subscore.labels = {};
            obj.subscore.pcavars = {};
            obj.subscore.weights = [];
            obj.subscore.colors{1,1} = ea_color_wes('all');
            obj.subscore.colors{1,2} = flip(ea_color_wes('all'));
            obj.subscore.vis.showposamount = repmat([25,25],10,1); %total of 10 subscores - will delete when we know the total number of subscores
            obj.subscore.vis.shownegamount = repmat([25,25],10,1);
            obj.subscore.vis.pos_shown = repmat([25,25],10,1);
            obj.subscore.vis.neg_shown = repmat([25,25],10,1);
            obj.subscore.negvisible = zeros(10,1);
            obj.subscore.posvisible = ones(10,1);
            obj.subscore.splitbysubscore = 0;
            obj.subscore.special_case = 0;
            obj.covarlabels={};
            obj.fiberdrawn = struct();
            obj.spotdrawn = struct();
            obj.networkdrawn = struct();
            obj.vizmode = 'Regions';
            obj.model = 'Smoothed';
            obj.smooth_fp = 0;
            obj.normalize_fp = 0;
            obj.mask_vta_fp = 0;
            obj.cvmask = 'Gray Matter';
            obj.fileformatversion=1.2; % new current version with settings to harmonize stats.
            datapath = GetFullPath(datapath);
            U = load(datapath, '-mat');
            if isfield(U, 'M') % Lead Group analysis path loaded
                obj.M = U.M;
                obj.leadgroup = datapath;

                testID = obj.M.guid;
                ea_mkdir([fileparts(obj.leadgroup),filesep,'UnifiedMappingExplorer',filesep]);
                id = 1;
                while exist([fileparts(obj.leadgroup),filesep,'UnifiedMappingExplorer',filesep,testID,'.explorer'],'file')
                    testID = [obj.M.guid, '_', num2str(id)];
                    id = id + 1;
                end
                obj.ID = testID;
                obj.resultfig = resultfig;

                if isfield(obj.M,'pseudoM')
                    obj.allpatients = obj.M.ROI.list;
                    obj.patientselection = 1:length(obj.M.ROI.list);
                    obj.M.root = [fileparts(datapath),filesep];
                    obj.M.patient.list = cell(size(obj.M.ROI.list,1), 1);
                    for i = 1:size(obj.M.ROI.list,1)
                        obj.M.patient.list{i,1} = obj.M.ROI.list{i,1};
                    end
                    obj.M.patient.group=obj.M.ROI.group; % copies
                else
                    datasetFolder = regexp(obj.leadgroup, ['(.*)(?=\', filesep, 'derivatives\', filesep, 'leadgroup)'], 'match', 'once');
                    for i = 1:size(obj.M.patient.list,1)
                        patient_tag = regexp(obj.M.patient.list{i}, '[^\\/]+$', 'match', 'once');
                        obj.M.patient.list{i} = fullfile(datasetFolder, 'derivatives', 'leaddbs', patient_tag);
                    end
                    obj.allpatients = obj.M.patient.list;
                    obj.patientselection = obj.M.ui.listselect;
                end
                obj.responsevar = obj.M.clinical.vars{1};
                obj.responsevarlabel = obj.M.clinical.labels{1};
                
            elseif  isfield(U, 'explorer')  % Saved explorer class loaded
                props = properties(U.explorer);
                for p =  1:length(props) %copy all public properties
                    if ~(strcmp(props{p}, 'analysispath') && ~isempty(obj.analysispath) ...
                            || strcmp(props{p}, 'ID') && ~isempty(obj.ID))
                        obj.(props{p}) = U.explorer.(props{p});
                    end
                end
                clear U
            else
                ea_error('You have opened a file of unknown type.')
                return
            end

            obj.repair_loaded_explorer;
            obj.compat_statmetric; % check and resolve for old statmetric code (which used to be integers)

            addlistener(obj,'activateby','PostSet',@activatebychange);

            % added a check here otherwise errors out for files w/o vatmodels
            if ~isfield(obj.M,'pseudoM')
                if ~isempty(obj.M.vatmodel) && contains(obj.M.vatmodel, 'OSS-DBS (Butenko 2020)')
                    obj.statmetric = 3;
                end
            end
        end

        function compat_statmetric(obj)
        
            if isempty(obj.statmetric) || ...
                    (isnumeric(obj.statmetric) && ~isscalar(obj.statmetric)) || ...
                    iscell(obj.statmetric)
                return
            end

            if ~ischar(obj.statmetric) % old language used:
                switch obj.statmetric % 3 was never used
                    case 1
                        obj.statmetric='Two-Sample T-Tests / VTAs (Baldermann 2019) / PAM (OSS-DBS)';
                    case 2
                        obj.statmetric='Correlations / E-fields (Irmen 2020)';
                    case 4
                        obj.statmetric='Proportion Test (Chi-Square) / VTAs (binary vars)';
                    case 5
                        obj.statmetric='Binomial Tests / VTAs (binary vars)';
                    case 6
                        obj.statmetric='Reverse T-Tests / E-Fields (binary vars)';
                    case 7
                        obj.statmetric='Plain Connections';
                    case 8
                        obj.statmetric='Odds Ratios / EF-Sigmoid (Jergas 2023)';
                    case 9
                        obj.statmetric='Weighted Linear Regression / EF-Sigmoid (Dembek 2023)';
                end
            end
            ea_unifiedmapping_compat_statmetrics2statsettings(obj);
        end


        function calculate(obj)
            % check that this has not been calculated before:
            %first store the rois for automatic calculations
            if ~isfield(obj.results,'roi') && ~(obj.calcsettings.connectivity_type == 2)
                if isfield(obj.M,'pseudoM')
                    vatlist = obj.M.ROI.list;
                else
                    vatlist = ea_sweetspotmapping_getvats(obj);
                end
                for vat=1:length(vatlist)
                    for side = 1:size(vatlist,2)
                        vta_nii = ea_load_nii(vatlist{vat,side});
                        obj.results.roi{vat,side} = vta_nii;
                    end
                end
            end
            if obj.calcsettings.selectedTool == 1 % sweetspotmapping
                  
                % in case of the sweetspot explorer, calculate rather means to
                % gather all E-Fields. To keep consistency of the logic with
                % discfiberexplorer and networkmappingexplorer, we will keep
                % the same name (calculate) for the function, nonetheless.

                % check that results aren't already there
                if isfield(obj.M,'pseudoM')
                    vatlist = obj.M.ROI.list;
                else
                    vatlist = ea_sweetspotmapping_getvats(obj);
                end


                [AllX,space] = ea_unifiedmapping_exportefieldmap(vatlist,obj);
                % Apply the calculation threshold while preserving
                % subthreshold observations as explicit zeros.
                for i = 1:numel(AllX)
                    if ~isempty(AllX{i})
                        AllX{i}(AllX{i} < obj.calcsettings.calcthreshold) = 0;
                    end
                end

                obj.results.sweetspotmapping.efield = AllX;
                
                obj.results.sweetspotmapping.space = space;

                if ~isfield(obj.M,'pseudoM')
                    try
                    % get active coordinates, as well
                    for pt=1:length(obj.M.patient.list)
                        for side=1:2
                            obj.results.sweetspotmapping.activecnt{side}(pt,:)=...
                                mean(obj.M.elstruct(pt).coords_mm{side}(find(obj.M.S(pt).activecontacts{side}),:),1); %#ok<FNDSB> % find is necessary here
                        end
                    end
                    for side=1:2
                        obj.results.sweetspotmapping.activecnt{side}=[obj.results.sweetspotmapping.activecnt{side};ea_flip_lr_nonlinear(obj.results.sweetspotmapping.activecnt{side})];
                    end
                    end
                end
            elseif obj.calcsettings.selectedTool == 2 % Fiber Filtering
                % if multi_pathways = 1, assemble cfile from multiple
                % pathway.dat files in dMRI_MultiTract/Connectome_name/
                % stores the result in the LeadGroup folder
                % also merges fiberActivation_.._.mat and stores them in
                % stimulation folders
                if ~isempty(obj.results) % something has been calculated
                    if isfield(obj.results,'fiberfiltering')
                        fname = ea_unifiedmapping_conn2connid(obj.calcsettings.fibfilt_connectome);
                        if isfield(obj.results.fiberfiltering, fname)
                            connField = obj.results.fiberfiltering.(fname);
                
                            % now check the specific subfields safely
                            if (isfield(connField,'PAM_Ttest') && obj.calcsettings.connectivity_type==2) || ...
                               (isfield(connField,'efield_mean') && obj.calcsettings.connectivity_type==1)
                
                                answ = questdlg('This has already been calculated. Are you sure you want to re-calculate everything?', ...
                                                'Recalculate Results','No','Yes','No');
                                if ~strcmp(answ,'Yes')
                                    return
                                end
                            end
                        end
                    end
                end
                if obj.calcsettings.multi_pathways == 1
                    [cfile, obj.map_list, obj.pathway_list] = ea_unifiedmapping_mergePathways(obj);
                else
                    cfile = [ea_getconnectomebase('dMRI'), obj.calcsettings.fibfilt_connectome, filesep, 'data.mat'];
                end
                %check if files exist
                FilesExist = check_stimvols(obj);

                if isfield(obj.M,'pseudoM') % failsave - this should not be necessary but still making sure things are set correctly for the pseudoM case.
                    obj.calcsettings.connectivity_type=1;
                    obj.calcsettings.calcmethod=1;
                end

                switch obj.calcsettings.connectivity_type
                    case 2    % if PAM, then just extracts activation states from fiberActivation.mat
                        fprintf("Calculating using the PAM method. Using dMRI connectome: %s",obj.calcsettings.fibfilt_connectome);
                        if all(FilesExist)
                            calculate_on_pam(obj,cfile)
                        end
                    otherwise     % check fiber recruitment via intersection with VTA
                        if obj.calcsettings.calcmethod == 1 %'E-field/Voxel Based Method'
                            fprintf("Calculating using the traditional E-field based method. Using dMRI connectome: %s",obj.calcsettings.fibfilt_connectome);
                            if all(FilesExist)
                                calculate_on_efield(obj,cfile)
                            end
                        elseif obj.calcsettings.calcmethod == 2 %'Fiber Based Method'
                            % check whether to use new (calc_on_fibers) or old method:
                            fprintf("Calculating using the Fiber based method. Using dMRI connectome: %s",obj.calcsettings.fibfilt_connectome);
                            if all(FilesExist)%why can this not be psuedo M?
                                calculate_on_fibers(obj,cfile)
                            end
                        end


                end

            elseif obj.calcsettings.selectedTool == 3 % network mapping
                
                if ~isempty(obj.results) % something has been calculated

                    if isfield(obj.results,'networkmapping')
                        if isfield(obj.results.networkmapping,ea_unifiedmapping_conn2connid(obj.calcsettings.netmap_connectome))
                            answ=questdlg('This has already been calculated. Are you sure you want to re-calculate everything?','Recalculate Results','No','Yes','No');
                            if ~strcmp(answ,'Yes')
                                return
                            end
                        end
                    end
                end
                
                if isfield(obj.M,'pseudoM')
                    vatlist = obj.M.ROI.list;
                else
                    %TODO:I have removed this from the networkmapping explorer folder and added it to the unified mapping explorer. Please adjust based on the future of the tool. I refrained from making a copy since the name of this script makes sense and would be redundant to change the name
                    vatlist = ea_unified_nm_getvats(obj); 
                end
                %TODO:I have removed this from the networkmapping explorer folder and added it to the unified mapping explorer. Please adjust based on the future of the tool. I refrained from making a copy since the name of this script makes sense and would be redundant to change the name
                [AllX, AllXVTAMasked] = ea_unified_nm_calcvals(vatlist, obj.calcsettings.netmap_connectome, obj.mask_vta_fp);

                obj.results.networkmapping.(ea_unifiedmapping_conn2connid(obj.calcsettings.netmap_connectome)).connval = AllX;
                if ~isempty(AllXVTAMasked)
                    obj.results.networkmapping.(ea_unifiedmapping_conn2connid(obj.calcsettings.netmap_connectome)).connval_vtamasked = AllXVTAMasked;
                elseif obj.mask_vta_fp
                    obj.mask_vta_fp = 0;
                end

                % Functional connectome, add spacedef to results
                if contains(obj.calcsettings.netmap_connectome, ' > ')
                    connid = (ea_unifiedmapping_conn2connid(obj.calcsettings.netmap_connectome));
                    connName = regexprep(obj.calcsettings.netmap_connectome, ' > .*$', '');
                    load([ea_getconnectomebase('fmri'), connName, filesep, 'dataset_volsurf.mat'], 'vol');
                    obj.results.networkmapping.(connid).space = vol.space;
                    if isfield(vol,'cifti')
                        obj.results.networkmapping.(connid).space.cifti=vol.cifti; % add cifti space as well (hidden in nii space)
                        obj.results.networkmapping.(connid).space.outidx=vol.outidx;
                        obj.results.networkmapping.(connid).space.inidx=vol.inidx;
                    end
                end
            end
            obj.hasResults = true;
        end




        function FilesExist = check_stimvols(obj)
            switch obj.calcsettings.connectivity_type
                case 2
                    [~,FilesExist] = ea_unifiedmapping_getpams(obj);
                    if ~all(FilesExist)
                        answ=questdlg('It seems like PAM has not been (completely) run. We can initiate the process now, but this will take some time. Proceed?','PAM not run','yes','no','yes');
                        switch answ
                            case 'yes'
                                options=ea_defaultoptions;
                                options.prefs.machine.vatsettings.butenko_calcPAM=1;
                                options.prefs.machine.vatsettings.butenko_calcVAT=0;
                                options.prefs.machine.vatsettings.butenko_connectome=obj.calcsettings.fibfilt_connectome;
                                options.groupdir=fileparts(obj.leadgroup);
                                obj.M.vatmodel='OSS-DBS (Butenko 2020)';
                                if isfield(obj.M.ui, 'stimSetMode') && obj.M.ui.stimSetMode
                                    options.stimSetMode = 1;
                                else
                                    options.stimSetMode = 0;
                                end
                                filesToCalc = find(sum(FilesExist(1:length(obj.M.patient.list),:),2)<2)';
                                calc_biophysical(obj,options,filesToCalc);
                            case 'no'
                                return
                        end

                    end
                    %recheck
                    [~,FilesExist] = ea_unifiedmapping_getpams(obj);
                otherwise
                    if obj.calcsettings.calcmethod == 2 %Fiber based method
                        [~,FilesExist] = ea_unifiedmapping_getlattice(obj);
                    else
                        if isfield(obj.M,'pseudoM')
                            for entry=1:length(obj.M.ROI.list)
                                FilesExist(entry)=exist(obj.M.ROI.list{entry},'file');
                            end
                        else
                            [~,FilesExist] = ea_unifiedmapping_getvats(obj);
                        end
                    end
                    while ~all(FilesExist(:))
                        answ=questdlg('It seems like not all stimulation volumes have been calculated. We can initiate the process now, but this will take some time. Proceed?','Stimvolumes not calculated','yes','no','yes');
                        switch answ
                            case 'yes'
                                if obj.calcsettings.calcmethod == 2 %Fiber based method
                                    obj.M.vatmodel='OSS-DBS (Butenko 2020)';
                                    switch obj.calcsettings.calcspace
                                        case 0
                                            space = 'native';
                                        case 1
                                            space = 'MNI';
                                    end
                                end
                                options=ea_defaultoptions;
                                options.prefs.machine.vatsettings.butenko_calcPAM=0;
                                options.prefs.machine.vatsettings.butenko_calcVAT=1;
                                options.groupdir=fileparts(obj.leadgroup);
                                if isfield(obj.M.ui, 'stimSetMode') && obj.M.ui.stimSetMode
                                    options.stimSetMode = 1;
                                else
                                    options.stimSetMode = 0;
                                end
                                filesToCalc = find(sum(FilesExist(1:length(obj.M.patient.list),:),2)<2)';
                                calc_biophysical(obj,options,filesToCalc);
                                [~,FilesExist] = ea_unifiedmapping_getvats(obj);
                            case 'no'
                                return
                        end
                    end
                    %recheck
                    
            end
            return
        end
        function calculate_on_pam(obj,cfile)
            connid = (ea_unifiedmapping_conn2connid(obj.calcsettings.fibfilt_connectome));
            [pamlist,~] = ea_unifiedmapping_getpams(obj);
            [fibsvalBin, fibsvalprob,~, ~, ~, fibcell_pam, connFiberInd, totalFibers] = ea_unifiedmapping_calcvals_pam_prob(pamlist, obj, cfile);

            obj.results.fiberfiltering.(connid).('PAM_probA').fibsval = fibsvalprob;
            obj.results.fiberfiltering.(connid).('PAM_Ttest').fibsval = fibsvalBin;
            obj.results.fiberfiltering.(connid).connFiberInd_PAM = connFiberInd;
            obj.results.fiberfiltering.(connid).totalFibers = totalFibers; % total number of fibers in the connectome to work with global indices
            obj.results.fiberfiltering.(connid).('pam_fibers').fibcell= fibcell_pam;
           
            % temp. duplicate fibcell, will be fixed in the new explorer
            obj.results.fiberfiltering.(ea_unifiedmapping_conn2connid(obj.calcsettings.fibfilt_connectome)).fibcell = obj.results.fiberfiltering.(ea_unifiedmapping_conn2connid(obj.calcsettings.fibfilt_connectome)).('pam_fibers').fibcell;
            %add a provision for the results 
            obj.results.fiberfiltering.(ea_unifiedmapping_conn2connid(obj.calcsettings.fibfilt_connectome)).calculationMethod = obj.calcsettings.calcmethod;
        end
        function calculate_on_efield(obj,cfile)
           connid = ea_unifiedmapping_conn2connid(obj.calcsettings.fibfilt_connectome);
           
            if isfield(obj.M,'pseudoM')
                vatlist=obj.M.ROI.list;
                [obj.customRoi.isbinary,obj.customRoi.minmax]=ea_unifiedmapping_checkcustomNii(vatlist);
                if obj.customRoi.isbinary
                    obj.statsettings.stimulationmodel='VTA';
                end
            else
                [vatlist,~] = ea_unifiedmapping_getvats(obj);
            end
            [fibsvalBin, fibsvalSum, fibsvalMean, fibsvalPeak, fibsval5Peak, fibcell_efield,  connFiberInd, totalFibers] = ea_fiberfiltering_calcvals(vatlist, cfile, obj.calcsettings.calcthreshold);
            obj.results.fiberfiltering.(connid).('VAT_Ttest').fibsval = fibsvalBin;
            obj.results.fiberfiltering.(connid).connFiberInd_VAT = connFiberInd; % old fiberfiltering files do not have these data and will fail when using pathway atlases
            obj.results.fiberfiltering.(connid).totalFibers = totalFibers; % total number of fibers in the connectome to work with global indices
            % only for e-fields
            obj.results.fiberfiltering.(connid).('efield_sum').fibsval = fibsvalSum;
            obj.results.fiberfiltering.(connid).('efield_mean').fibsval = fibsvalMean;
            obj.results.fiberfiltering.(connid).('efield_peak').fibsval = fibsvalPeak;
            obj.results.fiberfiltering.(connid).('efield_5peak').fibsval = fibsval5Peak;
            obj.results.fiberfiltering.(connid).('plainconn').fibsval = fibsvalBin;
            obj.results.fiberfiltering.(connid).('efield_fibers').fibcell= fibcell_efield;
            % temp. duplicate fibcell, will be fixed in the new explorer
            obj.results.fiberfiltering.(connid).fibcell = obj.results.fiberfiltering.(ea_unifiedmapping_conn2connid(obj.calcsettings.fibfilt_connectome)).('efield_fibers').fibcell;
            %add a provision for results
            obj.results.fiberfiltering.(connid).calculationMethod = obj.calcsettings.calcmethod;

        end

        function calculate_on_fibers(obj,cfile)
            connid = (ea_unifiedmapping_conn2connid(obj.calcsettings.fibfilt_connectome));
            
          
            % OSS-DBS E-field should be computed (not just warped!) in this space
            
            % get VAT list

            if isfield(obj.M,'pseudoM')
                vatlist = obj.M.ROI.list;
            else
                [vatlist,~] = ea_unifiedmapping_getlattice(obj);
            end
            
            
            % warp connectome to native space and compute E-field metrics
            ea_unified_get_Eproj(obj,vatlist)

            % define space again
            switch obj.calcsettings.calcspace
                case 1
                    space = 'MNI';
                case 0
                    space = 'native';
            end

            % load e-field projection metrics             
            [fibsvalBin_proj, fibsvalSum_proj, fibsvalMean_proj, fibsvalPeak_proj, fibsval5Peak_proj, fibcell_proj, connFiberInd_proj,fibsvalBin_magn, fibsvalSum_magn, fibsvalMean_magn, fibsvalPeak_magn, fibsval5Peak_magn, fibcell_magn, connFiberInd_magn, totalFibers] = ea_unifiedmapping_native_calcvals(vatlist, cfile, space, obj);

            obj.results.fiberfiltering.(connid).totalFibers = totalFibers; % total number of fibers in the connectome to work with global indices
            obj.results.fiberfiltering.(connid).('VAT_Ttest').fibsval = fibsvalBin_magn;
            obj.results.fiberfiltering.(connid).('efield_sum').fibsval = fibsvalSum_magn;
            obj.results.fiberfiltering.(connid).('efield_mean').fibsval = fibsvalMean_magn;
            obj.results.fiberfiltering.(connid).('efield_peak').fibsval = fibsvalPeak_magn;
            obj.results.fiberfiltering.(connid).('efield_5peak').fibsval = fibsval5Peak_magn;
            obj.results.fiberfiltering.(connid).('plainconn').fibsval = fibsvalBin_magn;
            obj.results.fiberfiltering.(connid).('efield_fibers').fibcell = fibcell_magn;
            obj.results.fiberfiltering.(connid).('efield_fibers').connFiberInd_VAT = connFiberInd_magn; % old fiberfiltering files do not have these data and will fail when using pathway atlases
            
            obj.results.fiberfiltering.(connid).('VAT_Ttest_proj').fibsval = fibsvalBin_proj;
            obj.results.fiberfiltering.(connid).('efield_proj_sum').fibsval = fibsvalSum_proj;
            obj.results.fiberfiltering.(connid).('efield_proj_mean').fibsval = fibsvalMean_proj;
            obj.results.fiberfiltering.(connid).('efield_proj_peak').fibsval = fibsvalPeak_proj;
            obj.results.fiberfiltering.(connid).('efield_proj_5peak').fibsval = fibsval5Peak_proj;
            obj.results.fiberfiltering.(connid).('plainconn_proj').fibsval = fibsvalBin_proj;
            obj.results.fiberfiltering.(connid).('efield_proj').fibcell = fibcell_proj;
            obj.results.fiberfiltering.(connid).('efield_proj').connFiberInd_VAT = connFiberInd_proj; % old fiberfiltering files do not have these data and will fail when using pathway atlases
            obj.results.fiberfiltering.(connid).calculationMethod = 'Fiber Based Method';
            if strcmp(obj.e_field_metric,'Magnitude')
                obj.results.fiberfiltering.(connid).fibcell = obj.results.fiberfiltering.(ea_unifiedmapping_conn2connid(obj.calcsettings.fibfilt_connectome)).('efield_fibers').fibcell;
                obj.results.fiberfiltering.(connid).connFiberInd_VAT = obj.results.fiberfiltering.(ea_unifiedmapping_conn2connid(obj.calcsettings.fibfilt_connectome)).('efield_fibers').connFiberInd_VAT;
            else
                obj.results.fiberfiltering.(connid).fibcell = obj.results.fiberfiltering.(ea_unifiedmapping_conn2connid(obj.calcsettings.fibfilt_connectome)).('efield_proj').fibcell;
                obj.results.fiberfiltering.(connid).connFiberInd_VAT = obj.results.fiberfiltering.(ea_unifiedmapping_conn2connid(obj.calcsettings.fibfilt_connectome)).('efield_proj').connFiberInd_VAT;
            end

        end

        function recalculate_fiberfiltering_threshold(obj)
            % Recomputes fiber connectivity with current obj.calcsettings.calcthreshold.
            % Call this after changing the E-field threshold so results and drawing use
            % the new threshold. Applies when Fiber Filtering is selected and E-field
            % (voxel-based) or Fiber-based method is used.
            % After calling, the GUI should refresh stats and redraw (e.g. obj.draw()).
            if obj.calcsettings.selectedTool ~= 2
                ea_cprintf('CmdWinWarnings', 'recalculate_fiberfiltering_threshold: Fiber Filtering is not selected. No action.\n');
                return;
            end
            connid = ea_unifiedmapping_conn2connid(obj.calcsettings.fibfilt_connectome);
            if ~isfield(obj.results,'fiberfiltering') || ~isfield(obj.results.fiberfiltering, connid)
                ea_cprintf('CmdWinWarnings', 'No prior fiberfiltering results for current connectome. Run Calculate first.\n');
                return;
            end
            if obj.calcsettings.multi_pathways == 1
                [cfile, obj.map_list, obj.pathway_list] = ea_unifiedmapping_mergePathways(obj);
            else
                cfile = [ea_getconnectomebase('dMRI'), obj.calcsettings.fibfilt_connectome, filesep, 'data.mat'];
            end
            FilesExist = check_stimvols(obj);
            if ~all(FilesExist(:))
                ea_cprintf('CmdWinWarnings', 'Not all stimulation volumes exist. Recalculation aborted.\n');
                return;
            end
            switch obj.calcsettings.calcmethod
                case 1
                    fprintf('Recalculating fiber connectivity with E-field threshold %g ...\n', obj.calcsettings.calcthreshold);
                    calculate_on_efield(obj, cfile);
                case 2
                    fprintf('Recalculating fiber connectivity (fiber-based) with E-field threshold %g ...\n', obj.calcsettings.calcthreshold);
                    calculate_on_fibers(obj, cfile);
                otherwise
                    ea_cprintf('CmdWinWarnings', 'recalculate_fiberfiltering_threshold: unsupported calcmethod. No action.\n');
            end
        end
        
        function  results = calc_biophysical(obj,options,filesToCalc)
            
            for pt = filesToCalc
                [options.root, options.patientname] = fileparts(obj.M.patient.list{pt});
                options.root = [options.root, filesep];
                options = ea_getptopts(fullfile(options.root, options.patientname), options);
                fprintf('\nProcessing %s...\n\n', options.patientname);

                if ~isfield(obj.M,'S')
                    ea_error(['Stimulation parameters for ', options.subj.subjId, ' are not set.']);
                end

                vfs = ea_regexpdir(ea_getearoot, 'ea_genvat_.*\.m$', 0);
                vfs = regexp(vfs, '(ea_genvat_.*)(?=\.m)', 'match', 'once');
                [vfnames,~,~] = cellfun(@(x) eval([x, '(''prompt'');']), vfs, 'Uni', 0);

                [~,ix]=ismember(obj.M.vatmodel,vfnames);
                try
                    ea_genvat=eval(['@',vfs{ix}]);
                catch
                    keyboard
                end
                if ~isfield(options.subj, 'norm')
                    ea_cprintf('CmdWinWarnings', 'Running in Miniset mode: %s...\n', options.subj.subjId);
                    volumespresent=0;
                elseif isempty(dir([options.subj.norm.transform.inverseBaseName, '*']))
                    ea_cprintf('CmdWinWarnings', 'Tranformation not found for %s...\n', options.subj.subjId);
                    volumespresent=0;
                else
                    volumespresent=1;
                end
                options.orignative=options.native; % backup
                options.native=~ea_getprefs('vatsettings.estimateInTemplate'); % see whether VTAs should be directly estimated in template space or not
                if options.native && ~volumespresent
                    ea_cprintf('CmdWinWarnings', 'Calculating VTA in template space since patient folder %s is incomplete.\n', options.subj.subjId);
                    options.native=0;
                end

                if options.native % Reload native space coordinates
                    coords = ea_load_reconstruction(options);
                else
                    coords = obj.M.elstruct(pt).coords_mm;
                end

                if strcmp(obj.M.vatmodel, 'OSS-DBS (Butenko 2020)')
                    if options.prefs.machine.vatsettings.butenko_calcAxonActivation
                        feval(ea_genvat,obj.M.S(pt),options);
                        ea_cprintf('CmdWinWarnings', 'OSS-DBS axon activation mode detect, skipping calc stats for %s!\n', options.patientname);
                        continue;
                    else
                        [vatCalcPassed, ~] = feval(ea_genvat,obj.M.S(pt),options);
                    end
                else
                    for side=1:2
                        try
                            [vtafv,vtavolume] = feval(ea_genvat,coords,obj.M.S(pt),side,options,['gs_',obj.M.guid]);
                            vatCalcPassed(side) = 1;
                        catch
                            vatCalcPassed(side) = 0;
                        end
                        if ~vatCalcPassed(side) 
                             ea_cprintf('CmdWinWarnings', 'VTA calculation failed for %s!\n', options.patientname);
                        end
                    end
                end

                options.native=options.orignative; % restore
            end
        end


        

        function Amps = getstimamp(obj)
            Amps=zeros(length(obj.M.patient.list),2);
            for pt=1:length(obj.M.patient.list)
                for side=1:2
                    thisamp=obj.M.stats(pt).ea_stats.stimulation.vat(side).amp;
                    thisamp(thisamp==0)=nan;
                    Amps(pt,side)=ea_nanmean(thisamp');
                end
            end
        end

        function VTAvolumes = getvtavolumes(obj)
            if ~isfield(obj.M.stats(1).ea_stats.stimulation.vat(1),'volume')
                VTAvolumes = obj.getstimamp;
                warning('No VTA volumes found. Using stimulation amplitudes instead. Re-run stats in Lead-group to obtain volumes.');
                return
            end
            VTAvolumes=zeros(length(obj.M.patient.list),2);
            for pt=1:length(obj.M.patient.list)
                for side=1:2
                    VTAvolumes(pt,side)=obj.M.stats(pt).ea_stats.stimulation.vat(side).volume;
                end
            end
        end

        function Efieldmags = getefieldmagnitudes(obj)
            if ~isfield(obj.M.stats(1).ea_stats.stimulation.efield(1),'volume')
                Efieldmags = obj.getstimamp;
                warning('No Efield magnitude sums found. Using stimulation amplitudes instead. Re-run stats in Lead-group to obtain values.');
                return
            end
            Efieldmags=zeros(length(obj.M.patient.list),2);
            for pt=1:length(obj.M.patient.list)
                for side=1:2
                    try
                        if isempty(obj.M.stats(pt).ea_stats.stimulation.efield(side).volume)
                            val=0;
                        else
                            val=obj.M.stats(pt).ea_stats.stimulation.efield(side).volume;
                        end
                        Efieldmags(pt,side)=val;
                    catch % could be efield(side) is not defined.
                        Efieldmags(pt,side)=0;
                    end
                end
            end
        end

        function refreshlg(obj)
            if ~exist(obj.leadgroup,'file')
                msgbox('Groupan alysis file has vanished. Please select file.');
                [fn,pth]=uigetfile();
                obj.leadgroup=fullfile(pth,fn);
            end
            U = load(obj.leadgroup);
            obj.M = U.M;
            obj.allpatients=obj.M.patient.list;
            obj.repair_loaded_explorer;
        end

        function repair_loaded_explorer(obj)
            % Fill fields that older/shared .explorer files may not carry.
            if isempty(obj.statmetric)
                obj.statmetric = nan;
            end
            if isempty(obj.statsettings) || ~isstruct(obj.statsettings)
                obj.statsettings = struct;
            end
            if ~isfield(obj.statsettings,'doVoxels') || isempty(obj.statsettings.doVoxels)
                obj.statsettings.doVoxels = 1;
            end
            if ~isfield(obj.statsettings,'doFibers') || isempty(obj.statsettings.doFibers)
                obj.statsettings.doFibers = 1;
            end
            if ~isfield(obj.statsettings,'outcometype') || isempty(obj.statsettings.outcometype)
                obj.statsettings.outcometype = 'gradual';
            end
            if ~isfield(obj.statsettings,'stimulationmodel') || isempty(obj.statsettings.stimulationmodel)
                obj.statsettings.stimulationmodel = 'Electric Field';
            end
            if ~isfield(obj.statsettings,'efieldmetric') || isempty(obj.statsettings.efieldmetric)
                obj.statsettings.efieldmetric = 'Sum';
            end
            if isfield(obj.statsettings, 'efieldthreshold') && ~isempty(obj.statsettings.efieldthreshold)
                legacyThreshold = obj.statsettings.efieldthreshold;
            else
                 legacyThreshold = 200;
            end
            if ~isfield(obj.statsettings, 'efieldthreshold_spot') || isempty(obj.statsettings.efieldthreshold_spot)
                obj.statsettings.efieldthreshold_spot = legacyThreshold;
            end

            if ~isfield(obj.statsettings,'efieldthreshold_tract') || isempty(obj.statsettings.efieldthreshold_tract)
                obj.statsettings.efieldthreshold_tract = legacyThreshold;
            end

            if ~isfield(obj.statsettings,'efieldthreshold_network') || isempty(obj.statsettings.efieldthreshold_network)
                obj.statsettings.efieldthreshold_network = 50;
            end
            
            if ~isfield(obj.statsettings,'nanthreshold')
                obj.statsettings.nanthreshold = [];
            end
            if ~isfield(obj.statsettings,'sweetspotresolution') || isempty(obj.statsettings.sweetspotresolution)
                obj.statsettings.sweetspotresolution = 0.5;
            end
            if ~isfield(obj.statsettings,'connthreshold') || isempty(obj.statsettings.connthreshold)
                obj.statsettings.connthreshold = 20;
            end
            if ~isfield(obj.statsettings,'statfamily') || isempty(obj.statsettings.statfamily)
                obj.statsettings.statfamily = 'Correlations';
            end
            if ~isfield(obj.statsettings,'stattest') || isempty(obj.statsettings.stattest)
                obj.statsettings.stattest = 'Spearman';
            end
            if ~isfield(obj.statsettings,'H0') || isempty(obj.statsettings.H0)
                obj.statsettings.H0 = 'Average';
            end

            if isempty(obj.calcsettings) || ~isstruct(obj.calcsettings)
                obj.calcsettings = struct;
            end
            if ~isfield(obj.calcsettings,'selectedTool') || isempty(obj.calcsettings.selectedTool)
                obj.calcsettings.selectedTool = 1;
            end
            if ~isfield(obj.calcsettings,'calcthreshold') || isempty(obj.calcsettings.calcthreshold)
                obj.calcsettings.calcthreshold = 50;
            end
            if ~isfield(obj.calcsettings,'switch_connectivity') || isempty(obj.calcsettings.switch_connectivity)
                obj.calcsettings.switch_connectivity = 1;
            end
            if ~isfield(obj.calcsettings,'connectivity_type') || isempty(obj.calcsettings.connectivity_type)
                obj.calcsettings.connectivity_type = 1;
            end
            if ~isfield(obj.calcsettings,'functionalresolution') || isempty(obj.calcsettings.functionalresolution)
                obj.calcsettings.functionalresolution = '2 mm';
            end
            if ~isfield(obj.calcsettings,'structuralresolution') || isempty(obj.calcsettings.structuralresolution)
                obj.calcsettings.structuralresolution = '2 mm';
            end
            if ~isfield(obj.calcsettings,'calcmethod') || isempty(obj.calcsettings.calcmethod)
                obj.calcsettings.calcmethod = 1;
            end
            if ~isfield(obj.calcsettings,'calcspace') || isempty(obj.calcsettings.calcspace)
                obj.calcsettings.calcspace = 1;
            end
            if ~isfield(obj.calcsettings,'netmap_connectome') || isempty(obj.calcsettings.netmap_connectome)
                obj.calcsettings.netmap_connectome = '';
            end
            if ~isfield(obj.calcsettings,'fibfilt_connectome') || isempty(obj.calcsettings.fibfilt_connectome)
                obj.calcsettings.fibfilt_connectome = '';
            end
            if ~isfield(obj.calcsettings,'multi_pathways') || isempty(obj.calcsettings.multi_pathways)
                obj.calcsettings.multi_pathways = 0;
            end

            if isempty(obj.subscore) || ~isstruct(obj.subscore)
                obj.subscore = struct;
            end
            if ~isfield(obj.subscore,'vars') || isempty(obj.subscore.vars)
                obj.subscore.vars = {};
            end
            if ~isfield(obj.subscore,'labels') || isempty(obj.subscore.labels)
                obj.subscore.labels = {};
            end
            if ~isfield(obj.subscore,'pcavars') || isempty(obj.subscore.pcavars)
                obj.subscore.pcavars = {};
            end
            if ~isfield(obj.subscore,'weights') || isempty(obj.subscore.weights)
                obj.subscore.weights = [];
            end
            if ~isfield(obj.subscore,'colors') || isempty(obj.subscore.colors)
                obj.subscore.colors{1,1} = ea_color_wes('all');
                obj.subscore.colors{1,2} = flip(ea_color_wes('all'));
            end
            if ~isfield(obj.subscore,'vis') || isempty(obj.subscore.vis)
                obj.subscore.vis = struct;
            end
            if ~isfield(obj.subscore.vis,'showposamount') || isempty(obj.subscore.vis.showposamount)
                obj.subscore.vis.showposamount = repmat([25,25],10,1);
            end
            if ~isfield(obj.subscore.vis,'shownegamount') || isempty(obj.subscore.vis.shownegamount)
                obj.subscore.vis.shownegamount = repmat([25,25],10,1);
            end
            if ~isfield(obj.subscore.vis,'pos_shown') || isempty(obj.subscore.vis.pos_shown)
                obj.subscore.vis.pos_shown = repmat([0,0],10,1);
            end
            if ~isfield(obj.subscore.vis,'neg_shown') || isempty(obj.subscore.vis.neg_shown)
                obj.subscore.vis.neg_shown = repmat([0,0],10,1);
            end
            if ~isfield(obj.subscore,'negvisible') || isempty(obj.subscore.negvisible)
                obj.subscore.negvisible = zeros(10,1);
            end
            if ~isfield(obj.subscore,'posvisible') || isempty(obj.subscore.posvisible)
                obj.subscore.posvisible = ones(10,1);
            end
            if ~isfield(obj.subscore,'splitbysubscore') || isempty(obj.subscore.splitbysubscore)
                obj.subscore.splitbysubscore = 0;
            end
            if ~isfield(obj.subscore,'special_case') || isempty(obj.subscore.special_case)
                obj.subscore.special_case = 0;
            end

            if isempty(obj.activated) || ~isstruct(obj.activated)
                obj.activated = struct;
            end
            if ~isfield(obj.activated,'sweetspotmapping')
                obj.activated.sweetspotmapping = 'Off';
            end
            if ~isfield(obj.activated,'fiberfiltering')
                obj.activated.fiberfiltering = 'Off';
            end
            if ~isfield(obj.activated,'networkmapping')
                obj.activated.networkmapping = 'Off';
            end

            if isempty(obj.NMviz) || ~isstruct(obj.NMviz)
                obj.NMviz = struct;
            end
            if ~isfield(obj.NMviz,'modelRH') || isempty(obj.NMviz.modelRH)
                obj.NMviz.modelRH = 1;
            end
            if ~isfield(obj.NMviz,'modelLH') || isempty(obj.NMviz.modelLH)
                obj.NMviz.modelLH = 1;
            end
            if isempty(obj.vizmode)
                obj.vizmode = 'Regions';
            end
            if isempty(obj.model)
                obj.model = 'Smoothed';
            end
            if ~iscell(obj.drawvals)
                obj.drawvals = {};
            end

            hasM = isstruct(obj.M);

            if isempty(obj.patientselection) && hasM
                if isfield(obj.M,'ui') && isfield(obj.M.ui,'listselect') && ~isempty(obj.M.ui.listselect)
                    obj.patientselection = obj.M.ui.listselect;
                elseif isfield(obj.M,'patient') && isfield(obj.M.patient,'list')
                    obj.patientselection = 1:numel(obj.M.patient.list);
                elseif isfield(obj.M,'ROI') && isfield(obj.M.ROI,'list')
                    obj.patientselection = 1:size(obj.M.ROI.list,1);
                end
            end

            if isempty(obj.allpatients) && hasM
                if isfield(obj.M,'patient') && isfield(obj.M.patient,'list')
                    obj.allpatients = obj.M.patient.list;
                elseif isfield(obj.M,'ROI') && isfield(obj.M.ROI,'list')
                    obj.allpatients = obj.M.ROI.list;
                end
            end

            if hasM && isfield(obj.M,'clinical') && isfield(obj.M.clinical,'labels') && ...
                    isfield(obj.M.clinical,'vars') && ~isempty(obj.M.clinical.labels)
                labels = obj.M.clinical.labels;
                if isempty(obj.responsevarlabel) || ...
                        ~(ischar(obj.responsevarlabel) || isstring(obj.responsevarlabel)) || ...
                        ~ismember(obj.responsevarlabel, labels)
                    obj.responsevarlabel = labels{1};
                end
                [isVar, ix] = ismember(obj.responsevarlabel, labels);
                if (isempty(obj.responsevar) || ~isVar) && isVar
                    obj.responsevar = obj.M.clinical.vars{ix};
                end
            end
        end

        function coh = getcohortregressor(obj)
            coh=ea_cohortregressor(obj.M.patient.group(obj.patientselection));
        end

        function [I, Ihat] = loocv(obj,silent)
            if ~exist('silent','var')
                silent=0;
            end
            rng(obj.rngseed);
            cvp = cvpartition(length(obj.patientselection), 'LeaveOut');
            [I, Ihat] = crossval(obj, cvp,[],0,silent);
        end

        function [I, Ihat] = lococv(obj,silent)
            if length(unique(obj.M.patient.group(obj.patientselection))) == 1
                ea_error(sprintf(['Only one cohort in the analysis.\n', ...
                    'Leave-One-Cohort-Out-validation not possible.']));
            end
            [I, Ihat] = crossval(obj, obj.M.patient.group(obj.patientselection),[],0,silent);
        end

        function [I, Ihat, val_struct] = kfoldcv(obj,silent)
            if ~exist('silent','var')
                silent=0;
            end
            I_iter = {};
            Ihat_iter = {};
            rng(obj.rngseed);
            iter = obj.kIter;
            if iter == 1
                cvp = cvpartition(length(obj.patientselection),'KFold',obj.kfold);
                [I,Ihat, val_struct] = crossval(obj,cvp,[],0,silent);
            else
                % plot some statistics over shuffles
                r_over_iter = zeros(iter,1);
                p_over_iter = zeros(iter,1);
                for i=1:iter
                    cvp = cvpartition(length(obj.patientselection), 'KFold', obj.kfold);
                    if ~silent
                        fprintf("Iterating fold set: %d",i)
                    end
                    [I_iter{i}, Ihat_iter{i},val_struct] = crossval(obj, cvp, [], 1,silent);
                    if ~silent
                        switch obj.multitractmode
                            case 'Split & Color By PCA'
                                disp("Fold Agreement is not evaluated for PCA")
                            otherwise
                                inx_nnan = find(isnan(I_iter{i}) ~= 1);
                                [r_over_iter(i),p_over_iter(i)]=ea_permcorr(I_iter{i}(inx_nnan),Ihat_iter{i}(inx_nnan),'spearman');
                        end
                    end
                end

                % check model agreement over shuffles using Sequential Rank Agreement
                % disabled for PCA
                switch obj.multitractmode
                    case 'Split & Color By PCA'
                        if ~silent
                            disp("Fold Agreement is not evaluated for PCA")
                        end
                    otherwise
                        if ~silent
                            r_Ihat = zeros(size(Ihat_iter,2));
                            for i = 1:size(r_Ihat,1)
                                for j = 1:size(r_Ihat,1)
                                    [r_Ihat(i,j),~]=ea_permcorr(Ihat_iter{i},Ihat_iter{j},'spearman');
                                end
                            end
                            % plot correlation matrix
                            figure('Name','Patient scores'' correlations','Color','w','NumberTitle','off')
                            imagesc(triu(r_Ihat));
                            title('Patient scores'' correlations over K-fold shuffles', 'FontSize', 16); % set title
                            colormap('bone');
                            cb = colorbar;
                            set(cb)

                            % plot r-vals over shuffles
                            p_above_05 = p_over_iter(find(p_over_iter>0.05),:);
                            p_above_01 = p_over_iter(find(p_over_iter>0.01),:);
                            h = figure('Name','Over-fold analysis','Color','w','NumberTitle','off');
                            g = ea_raincloud_plot(r_over_iter,'box_on',1);
                            a1=gca;
                            set(a1,'ytick',[])
                            a1.XLabel.String='Spearman''s R of model and clinical scores';

                            if min(r_over_iter) >= -0.9
                                r_lower_lim = min(r_over_iter) - 0.1;
                            else
                                r_lower_lim = -1.0;
                            end
                            if max(r_over_iter) <= 0.9
                                r_upper_lim = max(r_over_iter) + 0.1;
                            else
                                r_upper_lim = 1.0;
                            end

                            a1.XLim=([r_lower_lim r_upper_lim]);
                            text(0.25,0.9,['N(p>0.05) = ',sprintf('%d',length(p_above_05))],'FontWeight','bold','FontSize',14,'HorizontalAlignment','right','Units','normalized');
                            text(0.25,0.83,['N(p>0.01) = ',sprintf('%d',length(p_above_01))],'FontWeight','bold','FontSize',14,'HorizontalAlignment','right','Units','normalized');
                        end
                end

                % we should think about this part
                I_iter = cell2mat(I_iter);
                Ihat_iter = cell2mat(Ihat_iter);
                I = mean(I_iter,2,'omitnan');
                Ihat = mean(Ihat_iter,2,'omitnan');
            end
        end

        function [I, Ihat, val_struct] = lno(obj, Iperm, silent)
            if ~exist('silent','var')
                silent=0;
            end
            rng(obj.rngseed);
            cvp = cvpartition(length(obj.patientselection), 'resubstitution');
            if ~exist('Iperm', 'var') || isempty(Iperm)
                [I, Ihat, val_struct] = crossval(obj, cvp, [], [], silent);
            else
                [I, Ihat, val_struct] = crossval(obj, cvp, Iperm, [], silent);
            end
        end

        function [Improvement, Ihat, val_struct] = crossval(obj, cvp, Iperm, shuffle, silent)
            if ~exist('silent','var')
                silent=0;
            end
            if ~exist('shuffle','var') || isempty(shuffle)
                shuffle=0;
            end
            if isnumeric(cvp) % cvp is crossvalind
                cvIndices = cvp;
                cvID = unique(cvIndices);
                cvp = struct;
                cvp.NumTestSets = length(cvID);
                for i=1:cvp.NumTestSets
                    cvp.training{i} = cvIndices~=cvID(i);
                    cvp.test{i} = cvIndices==cvID(i);
                end
            end

            % Check if patients are selected in the custom training/test list
            if isempty(obj.customselection)
                patientsel = obj.patientselection;
            else
                patientsel = obj.customselection;
            end

            switch obj.multitractmode
                case 'Split & Color By PCA'
                    if ~exist('Iperm', 'var') || isempty(Iperm)
                        %Improvement = obj.subscore.vars;
                        for i=1:length(obj.subscore.vars)
                            Improvement{i} = obj.subscore.vars{i}(patientsel);
                        end
                    else
                        for i=1:length(obj.subscore.vars)
                            Improvement{i} = Iperm(patientsel,i);
                        end
                    end
                otherwise
                    if ~exist('Iperm', 'var') || isempty(Iperm)
                        Improvement = obj.responsevar(patientsel,:);
                    else
                        Improvement = Iperm(patientsel,:);
                    end
            end

            % Ihat is the estimate of improvements (not scaled to real improvements)
            if strcmp(obj.multitractmode,'Single Tract Analysis')
                Ihat = nan(length(patientsel),2);
                Ihat_train_global = nan(cvp.NumTestSets,length(patientsel),2);
            else
                Ihat = nan(length(patientsel),2,length(obj.subscore.vars));
                Ihat_train_global = nan(cvp.NumTestSets,length(patientsel),2,length(obj.subscore.vars));
            end
            if strcmp(obj.drawTool,'fiberfiltering')
                if obj.useExternalModel == true && ~strcmp(obj.ExternalModelFile, 'None')
                    S = load(obj.ExternalModelFile);
                    if ~strcmp(ea_unifiedmapping_method2methodid(obj),S.fibsvalType)
                        waitfor(msgbox('Change the Model Setup! See terminal'));
                        disp('The loaded model uses: ')
                        disp(S.fibsvalType)
                    end

                    fibsval = full(obj.results.fiberfiltering.(ea_unifiedmapping_conn2connid(obj.calcsettings.fibfilt_connectome)).(S.fibsvalType).fibsval);
                else
                    fibsval = full(obj.results.fiberfiltering.(ea_unifiedmapping_conn2connid(obj.calcsettings.fibfilt_connectome)).(ea_unifiedmapping_method2methodid(obj)).fibsval);
                end
            else
                fibsval = {};
            end


            % for nested LOO, store some statistics
            if obj.nestedLOO
                Abs_pred_error = zeros(cvp.NumTestSets, 1);
                Predicted_scores = zeros(length(patientsel), 1);
                Slope = zeros(cvp.NumTestSets, 1);
                Intercept = zeros(cvp.NumTestSets, 1);

            end

            
            for c=1:cvp.NumTestSets
                if cvp.NumTestSets ~= 1
                    if ~silent
                        fprintf(['\nIterating set: %0',num2str(numel(num2str(cvp.NumTestSets))),'d/%d\n'], c, cvp.NumTestSets);
                    end
                end

                if isobject(cvp)
                    training = cvp.training(c);
                    test = cvp.test(c);
                elseif isstruct(cvp)
                    training = cvp.training{c};
                    test = cvp.test{c};
                end

                % now do LOO within the training group
                if obj.nestedLOO
                    % use all patients, but outer loop left-out is always 0
                    if strcmp(obj.multitractmode,'Single Tract Analysis')
                        Ihat_inner = nan(length(patientsel),2);
                        Ihat_train_global_inner = nan(cvp.NumTestSets,length(patientsel),2);
                    else
                        Ihat_inner = nan(length(patientsel),2,length(obj.subscore.vars));
                        Ihat_train_global_inner = nan(cvp.NumTestSets,length(patientsel),2,length(obj.subscore.vars));
                    end
                    for test_i = 1:length(training)
                        training_inner = training;
                        training_inner(test_i) = 0;

                        % check if inner and outer left-out match
                        if all(training_inner == training)
                            continue
                        end

                        test_inner = logical(zeros(length(training), 1));
                        test_inner(test_i) = logical(training(test_i));

                        % updates Ihat_inner(test_inner)
                        if ~exist('Iperm', 'var') || isempty(Iperm)
                            [Ihat_inner, ~, ~] = ea_compute_unified_model(c,obj, fibsval, Ihat_inner, Ihat_train_global_inner, patientsel, training_inner, test_inner);
                        else
                            [Ihat_inner, ~, ~] = ea_compute_unified_model(c, obj, fibsval, Ihat_inner, Ihat_train_global_inner, patientsel, training_inner, test_inner,Iperm);
                        end
                    end

                    % fit the linear model based on inner loop fibscores
                    predictor=squeeze(ea_nanmean(Ihat_inner,2));
                    % iterating over all test_inner gives us training
                    mdl=fitglm(predictor(training),Improvement(training),lower(obj.predictionmodel));
                    Intercept(c) = mdl.Coefficients.Estimate(1);
                    Slope(c) = mdl.Coefficients.Estimate(2);
                end

                % now compute Ihat for the true 'test' left out
                % updates Ihat(test)
                if ~exist('Iperm', 'var') || isempty(Iperm)
                    [Ihat, Ihat_train_global, val_struct{c}] = ea_compute_unified_model(c,obj, fibsval, Ihat, Ihat_train_global, patientsel, training, test);
                else
                    [Ihat, Ihat_train_global, val_struct{c}] = ea_compute_unified_model(c,obj, fibsval, Ihat, Ihat_train_global, patientsel, training, test, Iperm);
                end

                % predict the improvement in the left-out patient (fold) of
                % the outer loop
                if obj.nestedLOO
                    predictor=squeeze(ea_nanmean(Ihat,2));
                    Ihat_voters_prediction = repmat(predict(mdl,predictor(test)),1,2);
                    %Abs_pred_error(c) = abs(Improvement(test) - Ihat_voters_prediction(test));
                    Predicted_scores(test) = Ihat_voters_prediction(1:end,1); % only one value here atm
                end

            end

            % check if binary variable and not permutation test
            if (~exist('Iperm', 'var') || isempty(Iperm)) && all(ismember(Improvement(:,1), [0,1])) && size(val_struct{c}.vals,1) == 1
                % average across sides. This might be wrong for capsular response.
                Ihat_av_sides = ea_nanmean(Ihat,2);
                if isobject(cvp)
                    % In-sample
                    AUC = ea_logit_regression(0 ,Ihat_av_sides, Improvement, 1:size(Improvement,1), 1:size(Improvement,1));
                elseif isstruct(cvp)
                    % actual training and test
                    Ihat_train_global_av_sides = ea_nanmean(Ihat_train_global,3); % in this case, dimens is (1, N, sides)
                    AUC = ea_logit_regression(Ihat_train_global_av_sides(training)', Ihat_av_sides, Improvement, training, test);
                end

            end

            if ~silent
                % plot patient score correlation matrix over folds
                if (~exist('shuffle', 'var')) || shuffle == 0 || isempty(shuffle)
                    if cvp.NumTestSets ~= 1 && (strcmp(obj.multitractmode,'Single Tract Analysis') || strcmp(obj.multitractmode,'Single Tract Analysis Button'))

                        % put training and test scores together
                        Ihat_combined = cell(1,cvp.NumTestSets);
                        %Ihat_combined = Ihat_train_global;
                        for c=1:cvp.NumTestSets
                            if isobject(cvp)
                                training = cvp.training(c);
                                test = cvp.test(c);
                            elseif isstruct(cvp)
                                training = cvp.training{c};
                                test = cvp.test{c};
                            end

                            Ihat_combined{c}(training,1) = Ihat_train_global(c,training,1)';
                            Ihat_combined{c}(test,1) = Ihat(test,1);
                        end

                        r_Ihat = zeros(size(Ihat_combined,2));

                        for i = 1:size(r_Ihat,1)
                            for j = 1:size(r_Ihat,1)
                                [r_Ihat(i,j),~]=ea_permcorr(Ihat_combined{i},Ihat_combined{j},'spearman');
                            end
                        end

                        figure('Name','Patient scores'' correlations','Color','w','NumberTitle','off')
                        imagesc(triu(r_Ihat)); % Display correlation matrix as an image
                        title('Patient scores'' correlations over folds', 'FontSize', 16); % set title
                        colormap('bone');
                        cb = colorbar;
                        % set(cb)

                    end
                end
            end
            if obj.nestedLOO
                % cvs = 'L-O-O-O';
                % h = ea_corrbox(Improvement,Predicted_dif_models,'permutation',{['Disc. Fiber prediction ',upper(cvs)],empiricallabel,fibscorelabel});
                LM_values_slope = [num2str(mean(Slope)) ' ' char(177) ' ' num2str(std(Slope))];
                LM_values_intercept = [num2str(mean(Intercept)) ' ' char(177) ' ' num2str(std(Intercept))];
                disp('Mean and STD for slopes and intercepts of LMs')
                disp(LM_values_slope)
                disp(LM_values_intercept)

                % visualize lms and CIs for 5-fold or less
                if cvp.NumTestSets < 6
                    groups_nested = zeros(length(Predicted_scores),1);
                    for group_idx = 1:cvp.NumTestSets
                        groups_nested(cvp.test(group_idx)) = group_idx;
                    end
                    side = 1;
                    plotName = 'Fitting of linear models for K-folds using nested LOO';
                    empiricallabel = 'Empirical score';
                    pred_label = 'Predicted score';
                    h=ea_corrbox(Improvement,Predicted_scores,'permutation',{['Disc. Fiber prediction ',plotName],empiricallabel,pred_label, plotName, LM_values_slope, LM_values_intercept},groups_nested);
                end
            end

            if obj.doactualprediction % repeat loops partly to fit to actual response variables:
                Ihat_voters_prediction=nan(size(Ihat));
                %add some warnings
                switch obj.multitractmode
                    case 'Single Tract Analysis'
                        if obj.useExternalModel && size(val_struct{c}.vals,1) > 1
                            ea_error("You can only use the Fit-to-Score feature with a Single Tract Analysis analysis model");
                        end
                    otherwise
                        if obj.useExternalModel
                            ea_error("You can only use the Fit-to-Score feature with Single Tract Analysis");
                        end
                end
                numVoters = size(val_struct{c}.vals,1);
                 for c=1:cvp.NumTestSets
                    if isobject(cvp)
                        training = cvp.training(c);
                        test = cvp.test(c);
                    elseif isstruct(cvp)
                        training = cvp.training{c};
                        test = cvp.test{c};
                    end
                    for voter=1:numVoters
                        switch obj.multitractmode
                            case 'Split & Color By Subscore'
                                if ~exist('Iperm', 'var') || isempty(Iperm)
                                    useI=obj.subscore.vars{voter}(patientsel);
                                else 
                                    % to be added by Nanditha
                                end 
                            case 'Split & Color By PCA'
                                if ~exist('Iperm', 'var') || isempty(Iperm)
                                    useI=obj.subscore.pcavars{voter}(patientsel);
                                else
                                    PCscores = ea_nanzscore(Iperm(patientsel, : ))*obj.subscore.pcacoeff;
                                    useI = PCscores(:, voter); 
                                end
                            otherwise
                                if ~exist('Iperm', 'var') || isempty(Iperm)
                                    useI=obj.responsevar(patientsel);
                                else 
                                    useI=Iperm(patientsel); 
                                end 
                        end

                        if size(useI,2)>1
                            ea_error('This has not been implemented for hemiscores.');
                        end

                        % these predictors are defined within the same fiberfiltering model
                        % of iteration 'c'
                        % do not get rid of the first dimension when it has size of 1
                        Ihat_train_global_av_sides = ea_nanmean(Ihat_train_global,3);
                        predictor_training = reshape(Ihat_train_global_av_sides, ...
                            size(Ihat_train_global_av_sides,1),...
                            size(Ihat_train_global_av_sides,2),...
                            size(Ihat_train_global_av_sides,4));

                        predictor_test = squeeze(ea_nanmean(Ihat,2));
                        %predictor=squeeze(ea_nanmean(Ihat_voters,2));

                        covariates=[];
                        for cv = 1:length(obj.covars)
                            covariates = [covariates,obj.covars{cv}(patientsel)];

                        end

                        if obj.useExternalModel == true %only use for single tract analysis
                            if ~strcmp(obj.multitractmode,'Single Tract Analysis')
                                ea_error("Sorry, you cannot use exported model and fit-to-scores for multi-tract model");
                            else
                                mdl = S.mdl;
                            end
                        else
                            if ~isempty(covariates)
                                mdl=fitglm([predictor_training(c,training,voter)',covariates(training,:)],useI(training),lower(obj.predictionmodel));
                            else
                                mdl=fitglm([predictor_training(c,training,voter)],useI(training),lower(obj.predictionmodel));
                            end
                        end
                        if size(useI,2) == 1 % global scores
                            if ~isempty(covariates)
                                Ihat_voters_prediction(test,:,voter)=repmat(predict(mdl,[predictor_test(test,voter),covariates(test,:)]),1,2); % fill both sides equally
                            else
                                Ihat_voters_prediction(test,:,voter)=repmat(predict(mdl,[predictor_test(test,voter)]),1,2); % fill both sides equally
                            end
                        elseif size(useI,2)==2 % bihemispheric scores
                            ea_error('Fitting to scores has not been implemented for bihemispheric scores.');
                        end
                    end
                end


                % quantify the prediction accuracy (if Train-Test)
                if cvp.NumTestSets == 1 && voter == 1 && size(obj.responsevar,2) == 1 && (~exist('Iperm', 'var') || isempty(Iperm))
                    side = 1;
                    SS_tot = var(useI(test)) * (length(useI(test)) - 1); % just a trick to use one line
                    SS_res = sum((Ihat_voters_prediction(test,side,1) - useI(test)).^2);
                    R2 = 1 - SS_res/SS_tot;
                    RMS = sqrt(mean((Ihat_voters_prediction(test,side,1) - useI(test)).^2));
                    MAD = median(abs(Ihat_voters_prediction(test,side,1) - useI(test)));
                    MAE = mean(abs(Ihat_voters_prediction(test,side,1) - useI(test)));

                    plotName = 'TRAIN-TEST';
                    R2_label = ['R2 = ', sprintf('%.3f',R2)];
                    RMS_label = ['RMS = ', sprintf('%.3f',RMS)];
                    MAD_label = ['MAD = ', sprintf('%.3f',MAD)];

                    empiricallabel = 'Empirical score';
                    pred_label = 'Predicted score';
                    h = ea_corrbox(useI(test),Ihat_voters_prediction(test,side,1),'permutation',{['Disc. Fiber prediction ',plotName],empiricallabel,pred_label, plotName, R2_label, RMS_label, MAD_label});
                    % h2 = ea_corrbox(-1*useI(test),Ihat_voters_prediction(test,side,1),'permutation',{['Disc. Fiber prediction ',plotName],empiricallabel,pred_label, plotName, R2_label, RMS_label, MAD_label});
                end

                Ihat=Ihat_voters_prediction; % replace with actual response variables.
            end

            switch obj.multitractmode
                case 'Split & Color By Subscore'
                    if ~obj.CleartuneOptim
                        % here we map back to the single response variable using a
                        % weightmatrix
                        if isempty(obj.customselection)
                            selected_pts = obj.patientselection;
                        else
                            selected_pts = obj.customselection;
                        end
                        weightmatrix=zeros(size(Ihat));
                        for voter=1:size(Ihat,3)
                            if ~isnan(obj.subscore.weights(voter)) % same weight for all subjects in that voter (slider was used)
                                weightmatrix(:,:,voter)=obj.subscore.weights(voter);
                            else % if the weight value is nan, this means we will need to derive a weight from the variable of choice
                                weightmatrix(:,:,voter)=repmat(ea_minmax(obj.subscore.weightvars{voter}(selected_pts)),1,size(weightmatrix,2)/size(obj.subscore.weightvars{voter}(selected_pts),2));
                                weightmatrix(:,:,voter)=weightmatrix(:,:,voter)./max(obj.subscore.weightvars{voter}(selected_pts)); % weight for unnormalized data across voters *
                                % *) e.g. in case one symptom - bradykinesia -
                                % has a max of 20, while a second - tremor -
                                % will have a max of 5, we want to equilize
                                % those. We use minmax() in the line above to
                                % get rid of negative values and use
                                % ./ea_nansum below to take the average.
                            end
                        end

                        for xx=1:size(Ihat,1) % make sure voter weights sum up to 1
                            for yy=1:size(Ihat,2)
                                %                     for xx=1:size(Ihat_voters,1) % make sure voter weights sum up to 1
                                %                         for yy=1:size(Ihat_voters,2)
                                weightmatrix(xx,yy,:)=weightmatrix(xx,yy,:)./ea_nansum(weightmatrix(xx,yy,:));
                            end
                        end
                        Ihat=ea_nansum(Ihat.*weightmatrix,3);
                    else
                        Ihat = Ihat(test,:,:);
                        Ihat = reshape(Ihat,2,length(obj.subscore.vars))';
                        Improvement = Improvement(test);
                        return;
                    end
                case 'Split & Color By PCA'

                    Ihat=squeeze(ea_nanmean(Ihat,2));
                    %Ihat_voters=squeeze(ea_nanmean(Ihat_voters,2)); % need to assume global scores here for now.

                    % map back to PCA:
                    for i=1:length(obj.subscore.vars)
                        selected_subscores{i} = obj.subscore.vars{i}(patientsel);
                    end
                    subvars=ea_nanzscore(cell2mat(selected_subscores));
                    if size(subvars,2) <= 2
                        ea_warndlg("You may not have enough subscores & this might result in errors. Please consider selecting more subscores.")
                    end

                    % [coeff,score,latent,tsquared,explained,mu]=pca(subvars,'Rows','complete');
                    % use saved weights to ensure consistency
                    coeff = obj.subscore.pcacoeff;

                    if ~silent
                        % show predictions for PC scores
                        if ~exist('Iperm', 'var') || isempty(Iperm) % avoid plotting for each permutation if using permutations!
                            for pcc=1:obj.numpcs
                                if obj.subscore.posvisible(pcc)==1 || obj.subscore.negvisible(pcc)==1 % don't try to plot if not showing any fibers for this PC
                                    ea_corrplot(obj.subscore.pcavars{pcc}(patientsel),Ihat(:,pcc), 'noperm', ...
                                        {['Disc. Fiber prediction for PC ',num2str(pcc)],'PC score (Empirical)','PC score (Predicted)'},...
                                        [], [], obj.subscore.pcacolors(pcc, :));
                                    % sum(obj.subscore.pcavars{pcc}(obj.patientselection) - score(:,pcc)) % quick check
                                end
                            end
                        end
                    end

                    % data is zscored, such as mu is 0 (+ some computer rounding error)
                    % then adding mean is not required
                    % also, we want to take scores of the chosen PCs ONLY,
                    % and multiply by coeff of these PCs (= how they map to
                    % the variables) to get estimated clinical scores
                    Ihatout = Ihat(:,1:obj.numpcs)*coeff(:,1:obj.numpcs)';
                    %Ihatout = Ihat*coeff(:,1:obj.numpcs)' + repmat(mu,size(score,1),1);
                    %Ihatout = Ihat_voters*coeff(:,1:obj.numpcs)' + repmat(mu,size(score,1),1);

                    Ihat = mat2cell(Ihatout, size(Ihatout,1), ones(1,length(obj.subscore.vars)));

                otherwise
                    Ihat=squeeze(Ihat);
                    %Ihat=squeeze(Ihat_voters);
            end
            if ~iscell(Ihat)
                if cvp.NumTestSets == 1
                    Ihat = Ihat(test,:);
                    Improvement = Improvement(test);
                end

                if size(obj.responsevar,2)==2 % hemiscores
                    Ihat = Ihat(:); % compare hemiscores (electrode wise)
                    Improvement = Improvement(:);
                else
                    Ihat = ea_nanmean(Ihat,2); % compare bodyscores (patient wise)
                end
            end

            % restore original view in case of live drawing
            if obj.cvlivevisualize
                obj.draw;
            end
        end


        function [Iperm, Ihat, R0, R1, pperm, Rp95, val_struct] = lnopb(obj, corrType, silent)
            if ~exist('corrType', 'var')
                corrType = 'Spearman';
            end
            if ~exist('silent','var')
                silent=0;
            end
            numPerm = obj.Nperm;

            if strcmp(obj.multitractmode,'Split & Color By PCA')
                Iperm = ea_shuffle(cell2mat(obj.subscore.vars'), numPerm, obj.patientselection, obj.rngseed);
                Iperm(2:numPerm+1,:,:) = Iperm;
                Iperm(1,:,:) = cell2mat(obj.subscore.vars');
                Ihat = cell(numPerm+1,1);

                R = zeros(numPerm+1, length(obj.subscore.vars));

                for perm=1:numPerm+1
                    if perm==1
                        if ~silent; fprintf('Calculating without permutation\n\n'); end
                        [~, Ihat{perm},thisval_struct] = lno(obj, [], silent);
                    else
                        if ~silent; fprintf('Calculating permutation: %d/%d\n\n', perm-1, numPerm); end
                        [~, Ihat{perm},thisval_struct] = lno(obj, squeeze(Iperm(perm,:,:)), silent);
                    end
                    val_struct{perm}=thisval_struct{1};
                    for subvar = 1:length(obj.subscore.vars)
                        R(perm,subvar) = corr(Iperm(perm, obj.patientselection, subvar)',...
                            Ihat{perm}{subvar},'type',corrType,'rows','pairwise');
                    end
                end

                R(isnan(R)) = 1e-5; % do not get rid of Nans

                % generate null distribution
                R1 = R(1,:);
                for subvar = 1:length(obj.subscore.vars)
                    R0(:,subvar) = sort(R(2:end,subvar), 'descend');
                    Rp95(subvar) = R0(round(0.05*numPerm),subvar);
                    pperm(subvar) = mean(abs(R0(:,subvar))>=abs(R1(subvar)));
                    if ~silent; fprintf(['Permuted p for ' obj.subscore.labels{subvar} ' = ' num2str(pperm(subvar)) '.\n']); end
                end

                % Return only selected I
                Iperm = Iperm(:,obj.patientselection,:);

            else % any mode except PCA
                Iperm = ea_shuffle(obj.responsevar, numPerm, obj.patientselection, obj.rngseed)';
                Iperm = [obj.responsevar, Iperm];
                Ihat = cell(numPerm+1, 1);

                R = zeros(numPerm+1, 1);

                for perm=1:numPerm+1
                    if perm==1
                        if ~silent; fprintf('Calculating without permutation\n\n'); end
                        [~, Ihat{perm},val_struct{perm}] = lno(obj, [], silent);
                    else
                        if ~silent; fprintf('Calculating permutation: %d/%d\n\n', perm-1, numPerm); end
                        [~, Ihat{perm},val_struct{perm}] = lno(obj, Iperm(:, perm), silent);
                    end

                    R(perm) = corr(Iperm(obj.patientselection,perm),Ihat{perm},'type',corrType,'rows','pairwise');
                end

                R(isnan(R)) = 1e-5;

                % generate null distribution
                R1 = R(1);
                R0 = sort((R(2:end)),'descend');
                Rp95 = R0(round(0.05*numPerm));
                pperm = mean(abs(R0)>=abs(R1));
                if ~silent; disp(['Permuted p = ',sprintf('%0.2f',pperm),'.']); end

                % Return only selected I
                Iperm = Iperm(obj.patientselection,:);
            end
        end

        function save(obj, saveas)
            obj.repair_loaded_explorer;

            % Create a temporary object with only the required fields
            
            % Get all properties of the object
            explorer = ea_unifiedmapping;
            Incprops = {'results','calcsettings','statsettings','leadgroup','ID','M', ...
                'subscore','responsevar','responsevarlabel','patientselection', ...
                'allpatients','activated','multitractmode','posvisible','negvisible', ...
                'showposamount','shownegamount','NMviz','vizmode','model'};
            for i = 1:length(Incprops)
               explorer.(Incprops{i}) = obj.(Incprops{i});
            end

            % Determine save path
            if nargin < 2 || isempty(saveas)
                % Original behaviour: auto-derive path from leadgroup
                if isempty(obj.analysispath)
                    [pth,~,~] = fileparts(obj.leadgroup);
                    obj.analysispath = [pth, filesep, 'UnifiedMappingExplorer', ...
                                        filesep, obj.ID, '.explorer'];
                    ea_mkdir([pth, filesep, 'UnifiedMappingExplorer']);
                end
                savepath = obj.analysispath;
            else
                % Save As: use the user-chosen path, and update analysispath
                savepath = saveas;
                obj.analysispath = saveas;
            end
        
            rf = obj.resultfig;   % stash fig handle before saving
            rd = obj.drawobject;  % stash drawing handle before saving
            try
                setappdata(rf, ['dt_', explorer.ID], rd);
            end
        
            save(savepath, 'explorer', '-v7.3');
            saveObjectToJson(obj);
            obj.resultfig = rf;
            obj.drawobject = rd;
          

            % 
            % %This is necessary to match the settings file
            % %only save results in this
            % if isempty(obj.analysispath)
            %     [pth,~,~] = fileparts(obj.leadgroup);
            %     obj.analysispath=[pth,filesep,'UnifiedMappingExplorer',filesep,obj.ID,'.explorer'];
            %     ea_mkdir([pth,filesep,'UnifiedMappingExplorer']);
            % end
            % rf=obj.resultfig; % need to stash fig handle for saving.
            % rd=obj.drawobject; % need to stash handle of drawing before saving.
            % try % could be figure is already closed.
            %     setappdata(rf,['dt_',explorer.ID],rd); % store handle of tract to figure.
            % end
            % 
            % save(obj.analysispath,'explorer','-v7.3');
            % saveObjectToJson(obj);
            % obj.resultfig=rf;
            % obj.drawobject=rd;
        end

        function saveObjectToJson(obj)
            % Convert object to a struct (including nested objects)
            
            voxtractsettings = objectToStruct(obj);

            % % % force setselection to be stored as a cell array otherwise json will get
            % % % it wrong
            % % voxtractsettings.setselections = struct('type','celllogical', ...
            % %     'data',{obj.setselections});            
            % Convert struct to JSON
            jsonStr = jsonencode(voxtractsettings, 'PrettyPrint', true);
            %define filepaths
            if isempty(obj.analysispath)
                [DBSMappingfolder,~,~] = fileparts(obj.leadgroup);
            else
                [DBSMappingfolder,~,~] = fileparts(obj.analysispath);
            end
            conn_val = 'default';
            switch obj.drawTool
                case 'sweetspotmapping'
                conn_val = 'default';
                case 'fiberfiltering'
                conn_val = ea_unifiedmapping_conn2connid(obj.calcsettings.fibfilt_connectome);
                case 'networkmapping'
                conn_val = ea_unifiedmapping_conn2connid(obj.calcsettings.netmap_connectome);
            end
        
            if ~isfolder(DBSMappingfolder)
                ea_mkdir(DBSMappingfolder)
            end

            % check for custom save path
            if isprop(obj, 'AdditionalSettingsSavePath') && ...
                    ~isempty(obj.AdditionalSettingsSavePath)
                customPath = obj.AdditionalSettingsSavePath;
                % enforce .json extension
                [folder, name, ext] = fileparts(customPath);
                name = ['Settings-', name];
                if isempty(ext)
                    ext = '.json';
                end
            
                if ~strcmpi(ext, '.json')
                    error('AdditionalSettingsSavePath must point to a .json file');
                end

                % if no folder provided, use DBSMappingfolder
                if isempty(folder)
                    folder = DBSMappingfolder;
                end
            
                jsonPath = fullfile(folder, [name ext]);
            
            else
                % default behavior
                jsonPath = fullfile(DBSMappingfolder, ...
                    ['Settings-', obj.ID, '.json']);
            end

            % jsonPath=[DBSMappingfolder,filesep,'Settings-',obj.ID,'_conn-',conn_val,'.json'];
            % 
            % Write JSON to a file
            fileID = fopen(jsonPath, 'w');
            if fileID == -1
                error('Cannot open file for writing.');
            end
            fprintf(fileID, '%s', jsonStr);
            fclose(fileID);
        end

        function s = objectToStruct(obj)
            % Convert an object to a struct, handling nested objects
            if nargin < 2
                ignoreList = {'results','resultfig','drawobject','M'}; %M should be present in the explorer mat file. This is because there are some complicated structures in M files that are not well translated in struct (for json encoding). % Default: Do not ignore any properties unless specified
            end

            if isobject(obj)
                props = properties(obj);
                s = struct();
                for i = 1:length(props)
                    propName = props{i};
                    propValue = obj.(props{i});
                    if ismember(propName, ignoreList) %skip some of the properties.
                        continue;
                    end
                    if isa(propValue, 'matlab.ui.Figure') || isa(propValue, 'handle') %also skip handles
                        continue;
                    end
                    if isobject(propValue)  % Recursively convert nested objects
                        s.(props{i}) = objectToStruct(propValue);
                    else
                        s.(props{i}) = propValue;
                    end
                end
            else
                s = obj; % If it's not an object, return as is (handles arrays, numbers, strings)
            end
        end

        function draw(obj)
            obj.repair_loaded_explorer;
            
            if ~isfield(obj.activated,'sweetspotmapping')
                obj.activated.sweetspotmapping='Off';
            end
            if ~isfield(obj.activated,'fiberfiltering')
                obj.activated.fiberfiltering='Off';
            end
            if ~isfield(obj.activated,'networkmapping')
                obj.activated.networkmapping='Off';
            end

            % delete prior spots:
            if isfield(obj.drawobject,'sweetspotmapping')
                if ~isempty(obj.drawobject.sweetspotmapping)
                    for s=1:numel(obj.drawobject.sweetspotmapping)
                        for ins=1:numel(obj.drawobject.sweetspotmapping{s})
                            try delete(obj.drawobject.sweetspotmapping{s}{ins}.toggleH); end
                            try delete(obj.drawobject.sweetspotmapping{s}{ins}.patchH); end
                            try delete(obj.drawobject.sweetspotmapping{s}{ins}); end
                        end
                    end
                end
            end
            % plot new spots
            switch lower(obj.activated.sweetspotmapping)
                case 'on'
                    obj.drawTool='sweetspotmapping';
                    ea_unified_draw(obj);
            end
            % delete prior tracts
            if isfield(obj.drawobject,'fiberfiltering')
                if ~isempty(obj.drawobject.fiberfiltering)
                    for s=1:numel(obj.drawobject.fiberfiltering)
                        surfArray=obj.drawobject.fiberfiltering{s};
                        for ins=1:numel(surfArray)
                            try delete(surfArray(ins).toggleH); end
                            try delete(surfArray(ins).patchH); end
                            try delete(surfArray(ins)); end
                        end
                    end
                    obj.drawobject.fiberfiltering{s} = [];
                end
            end
            % plot new tracts
            switch lower(obj.activated.fiberfiltering)
                case 'on'
                    obj.drawTool='fiberfiltering';
                    ea_unified_draw(obj);
            end
            % delete prior nets
            if isfield(obj.drawobject,'networkmapping')
                if ~isempty(obj.drawobject.networkmapping)
                    for s=1:numel(obj.drawobject.networkmapping)
                        for ins=1:numel(obj.drawobject.networkmapping{s})
                            try delete(obj.drawobject.networkmapping{s}{ins}.toggleH); end
                            try delete(obj.drawobject.networkmapping{s}{ins}.patchH); end
                            try delete(obj.drawobject.networkmapping{s}{ins}); end
                        end
                    end
                end
            end
            % plot new nets
            switch lower(obj.activated.networkmapping)
                case 'on'
                    obj.drawTool='networkmapping';
                    ea_unified_draw(obj);
            end

   
        end
    end
    methods (Static)
        function obj = loadobj(obj)
            if isstruct(obj)
                saved = obj;
                obj = ea_unifiedmapping;
                fields = fieldnames(saved);
                for f = 1:numel(fields)
                    if isprop(obj, fields{f})
                        obj.(fields{f}) = saved.(fields{f});
                    end
                end
            end
            obj.repair_loaded_explorer;
        end

        function changeevent(~,event)
            update_trajectory(event.AffectedObject,event.Source.Name);
        end
    end
end



function activatebychange(~,event)
      % activate_tractset();
end
function calculateIntersection(obj)
for nroi = 1:length(obj.roiintersectdata)
    vat = ea_load_nii(obj.roiintersectdata{nroi}); %use only one, otherwise drawing doesn't make sense
    thresh = obj.roithresh;
    vatInd = find(abs(vat.img(:))>thresh);
    [xvox, yvox, zvox] = ind2sub(size(vat.img), vatInd);
    vatmm = ea_vox2mm([xvox, yvox, zvox], vat.mat);
    for side = 1:2
        for i=1:size(obj.drawobject,1)
            vals = {};
            valsPeak = {};
            connected = [];
            trimmedFiberInd = [];
            resultFibers = obj.fiberdrawn.fibcell{i,side};
            if isempty(resultFibers)
                continue
            end
            fibers=ea_fibcell2fibmat(resultFibers);
            filter = all(fibers(:,1:3)>=min(vatmm),2) & all(fibers(:,1:3)<=max(vatmm), 2);
            if ~any(filter)
                zeros_arr = zeros(size(obj.drawobject{i,side},1),1);
                normwts = mat2cell(zeros_arr,ones(size(obj.drawobject{i,side},1),1));
                [obj.drawobject{i,side}.FaceAlpha]=normwts{:};
                continue
            end
            trimmedFiber = fibers(filter,:);
            % Map mm connectome fibers into VAT voxel space
            [trimmedFiberInd, ~, trimmedFiberID] = unique(trimmedFiber(:,4), 'stable');
            fibVoxInd = splitapply(@(fib) {ea_mm2uniqueVoxInd(fib, vat)}, trimmedFiber(:,1:3), trimmedFiberID);
            % Remove outliers
            fibVoxInd(cellfun(@(x) any(isnan(x)), fibVoxInd)) = [];
            trimmedFiberInd(cellfun(@(x) any(isnan(x)), fibVoxInd)) = [];
            connected = cellfun(@(fib) any(ismember(fib, vatInd)), fibVoxInd);
            vals = cellfun(@(fib) vat.img(intersect(fib, vatInd)), fibVoxInd(connected), 'Uni', 0);
            valsPeak{1}(trimmedFiberInd(connected)) = cellfun(@mean, vals);
            wts = cell2mat(valsPeak)';
            if ~isempty(wts)
                if length(wts) ~= size(obj.drawobject{i,side},1)
                    diff = length(wts) - size(obj.drawobject{i,side},1);
                    if diff < 0
                        wts = [wts;zeros(abs(diff),1)];
                    end
                end
                normwts = normalize(ea_contrast(wts,10,0),'range');
                normwts =  mat2cell(normwts,ones(size(normwts,1),1));
                if ~isempty(normwts) && ~isempty(obj.drawobject{i,side})
                    try
                        [obj.drawobject{i,side}.FaceAlpha]=normwts{:};
                        disp(['Changed alpha of tract',num2str(i)]);
                        normwts = {};
                    end
                end
            else %if it is not connected then they should have zero alpha!!
                zeros_arr = zeros(size(obj.drawobject{i,side},1),1);
                normwts = mat2cell(zeros_arr,ones(size(obj.drawobject{i,side},1),1));
                [obj.drawobject{i,side}.FaceAlpha]=normwts{:};
            end
        end
    end
end
end



function check_and_update_visibility(obj, field, ~, condition, group)
if nargin < 5  % If group is not provided, operate on obj directly
    group = [];
end

if eval(sprintf('obj.%s%s && all(values %s)', field, group_access(group), condition))
    eval(sprintf('obj.%s%s = 0;', field, group_access(group)));
    fprintf('\n')
    warning('off', 'backtrace');
    warning('No %s values found, %s is set to 0 now!', condition_description(condition), field);
    warning('on', 'backtrace');
    fprintf('\n')
end
end

function access = group_access(group)
if isempty(group)
    access = '';
else
    access = sprintf('(group)');
end
end

function desc = condition_description(condition)
if strcmp(condition, '<0')
    desc = 'positive';
else
    desc = 'negative';
end
end
function fibers=ea_fibcell2fibmat(fibers)
[idx,~]=cellfun(@size,fibers);
fibers=cell2mat(fibers);
idxv=zeros(size(fibers,1),1);
lid=1; cnt=1;
for id=idx'

    idxv(lid:lid+id-1)=cnt;
    lid=lid+id;
    cnt=cnt+1;
end
fibers=[fibers,idxv];
end
