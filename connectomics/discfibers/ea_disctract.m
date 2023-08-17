classdef ea_disctract < handle
    % Discriminative fiber class to handle visualizations of discriminative fibers in lead dbs resultfig / 3D Matlab figures
    % A. Horn

    properties (SetObservable)
        M % content of lead group project
        resultfig % figure handle to plot results
        ID % name / ID of discriminative fibers object
        calcthreshold % initial hard threshold to impose on (absolute) nifti files only when calculating the data
        posvisible = 1 % pos tract visible
        negvisible = 0 % neg tract visible
        roivisible = 0 % show ROI (usually VTAs)
        connfibvisible = 0 % show all connected tracts in white
        showposamount = [25 25] % two entries for right and left
        shownegamount = [25 25] % two entries for right and left
        connthreshold = 20
        efieldthreshold = 200
        statmetric = 'Correlations / E-fields (Irmen 2020)' % the following stat metric are available:
        % ’Two-Sample T-Tests / VTAs (Baldermann 2019) / PAM (OSS-DBS)’
        % ’One-Sample Tests / VTAs (Baldermann 2019) / PAM (OSS-DBS)’
        % ‘Correlations / E-fields (Irmen 2020)’
        % ‘Proportion Test (Chi-Square) / VTAs (binary vars)’
        % ‘Binomial Tests / VTAs (binary vars)’
        % ‘Reverse T-Tests / E-Fields (binary vars)’
        % ‘Plain Connections’
        threshstrategy = 'Percentage Relative to Peak'; % can be 'Relative to Amount' or 'Fixed Amount'
        multi_pathways = 0 % if structural connectome is devided into pathways (multiple .mat in dMRI_MultiTract)
        map_list % list that contains global indices of the first fibers in each pathway (relevant when multi_pathways = 1)
        pathway_list % list that contains names of pathways (relevant when multi_pathways = 1
        connFiberInd % legacy
        connectivity_type = 1 % 1 - VAT, 2 - PAM
        switch_connectivity = 0 % flag if connectivity type was changed in the GUI
        nestedLOO = false       % if true, will conducted LOO in the training set
        use_adjacency = false   % if true, 'smoothes' the fiber model by addings scores of fibers next to the pivotal ones
        ADJ = false             % contains the large but sparse adjacency matrix
        adj_scaler = 0.025 % coefficient that determines the degree of smoothing. Atm, hardwired, but could adaptive
        corrtype = 'Spearman' % correlation strategy in case of statmetric ==  Correlations / E-fields (Irmen 2020). In case of one-sample tests used for 'T-Tests' vs 'Wicoxon Signed Rank Tests'.
        efieldmetric = 'Peak' % if statmetric == ‘Correlations / E-fields (Irmen 2020)’, efieldmetric can calculate sum, mean or peak along tracts
        poscolor = [0.9176,0.2000,0.1373] % positive main color
        negcolor = [0.2824,0.6157,0.9725] % negative main color
        multitractmode = 'Single Tract Analysis' % multi mode now coded by this value
        numpcs = 4; % standard value of how many PCs to compute in case of PCA mode
        doactualprediction = 0; % set up nested CVs to carry out actual predictions of response variables
        predictionmodel = 'Linear'; % type of glm used to fit fiber values to actual scores
        showsignificantonly = 0
        alphalevel = 0.05
        multcompstrategy = 'FDR'; % could be 'Bonferroni'
        subscore
        currentune
        results
        % Subfields:
        % results.(connectomename).fibcell: cell of all fibers connected, sorted by side
        % results.(connectomename).ttests.fibsval % connection status for each fiber to each VTA
        % results.(connectomename).spearman_sum.fibsval % connection weights for each fiber to each VTA
        % results.(connectomename).spearman_mean.fibsval % connection weights for each fiber to each VTA
        % results.(connectomename).spearman_peak.fibsval % connection weights for each fiber to each VTA
        % results.(connectomename).spearman_5peak.fibsval % connection weights for each fiber to each VTA
        cleartuneresults % copy of results for auto tuning functions
        cleartuneefields % efields used to calc results
        cleartuneinjected % status to report file has injected values
        CleartuneOptim = 0;
        cleartunevars
        activateby={}; % entry to use to show fiber activations
        cvlivevisualize = 0; % if set to 1 shows crossvalidation results during processing.
        basepredictionon = 'Mean of Scores';
        fiberdrawn % struct contains fibercell and vals drawn in the resultfig
        drawobject % actual streamtube handle
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
        connectome % redundancy protocol only, name of underlying connectome
        colorbar % colorbar information
        useExternalModel = false
        ExternalModelFile = 'None'
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
        % misc
        runwhite = 0; % flag to calculate connected tracts instead of stat tracts
        SigmoidTransform = 0;  % flag to transform E-field to Sigmoids
        twoSampleWeighted = 1;  % flag for two-sample weighted linear regression
    end

    properties (Access = private)
        switchedFromSpace=3 % if switching space, this will protocol where from
    end

    methods
        function obj=ea_disctract(analysispath) % class constructor
            if exist('analysispath', 'var') && ~isempty(analysispath)
                obj.analysispath = analysispath;
                [~, ID] = fileparts(obj.analysispath);
                obj.ID = ID;
            end
        end

        function initialize(obj,datapath,resultfig)

            datapath = GetFullPath(datapath);
            D = load(datapath, '-mat');
            if isfield(D, 'M') % Lead Group analysis path loaded
                obj.M = D.M;
                obj.leadgroup = datapath;

                testID = obj.M.guid;
                ea_mkdir([fileparts(obj.leadgroup),filesep,'fiberfiltering',filesep]);
                id = 1;
                while exist([fileparts(obj.leadgroup),filesep,'fiberfiltering',filesep,testID,'.fibfilt'],'file')
                    testID = [obj.M.guid, '_', num2str(id)];
                    id = id + 1;
                end
                obj.ID = testID;
                obj.resultfig = resultfig;

                if isfield(obj.M,'pseudoM')
                    obj.allpatients = obj.M.ROI.list;
                    obj.patientselection = 1:length(obj.M.ROI.list);
                    obj.M = ea_map_pseudoM(obj.M);
                    obj.M.root = [fileparts(datapath),filesep];
                    obj.M.patient.list = cell(size(obj.M.ROI.list,1), 1);
                    for i = 1:size(obj.M.ROI.list,1)
                        obj.M.patient.list{i,1} = obj.M.ROI.list{i,1};
                    end
                    obj.M.patient.group=obj.M.ROI.group; % copies
                else
                    obj.allpatients = obj.M.patient.list;
                    obj.patientselection = obj.M.ui.listselect;
                end

                obj.responsevarlabel = obj.M.clinical.labels{1};
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
            elseif  isfield(D, 'tractset')  % Saved tractset class loaded
                props = properties(D.tractset);
                for p =  1:length(props) %copy all public properties
                    if ~(strcmp(props{p}, 'analysispath') && ~isempty(obj.analysispath) ...
                            || strcmp(props{p}, 'ID') && ~isempty(obj.ID))
                        obj.(props{p}) = D.tractset.(props{p});
                    end
                end
                clear D
            else
                ea_error('You have opened a file of unknown type.')
                return
            end
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
        end
        function calculate(obj)

            switch obj.connectivity_type
                case 2 % PAM
                    obj.M.vatmodel = 'OSS-DBS (Butenko 2020)';
                otherwise % Stim. volumes
                    obj.M.vatmodel= 'SimBio/FieldTrip (see Horn 2017)';
            end

            % check that this has not been calculated before:
            if ~isempty(obj.results) % something has been calculated
                if isfield(obj.results,ea_conn2connid(obj.connectome))
                    answ=questdlg('This has already been calculated. Are you sure you want to re-calculate everything?','Recalculate Results','No','Yes','No');
                    if ~strcmp(answ,'Yes')
                        return
                    end
                end
            end

            % if multi_pathways = 1, assemble cfile from multiple
            % pathway.dat files in dMRI_MultiTract/Connectome_name/
            % stores the result in the LeadGroup folder
            % also merges fiberActivation_.._.mat and stores them in
            % stimulation folders
            if obj.multi_pathways == 1
                [cfile, obj.map_list, obj.pathway_list] = ea_discfibers_merge_pathways(obj);
            else
                cfile = [ea_getconnectomebase('dMRI'), obj.connectome, filesep, 'data.mat'];
            end

            % if multi_pathway = 1, use the adjacency defined for
            % merged_pathways.mat
            if obj.use_adjacency
                if obj.multi_pathways == 1
                    connectome_folder = [ea_getconnectomebase('dMRI_multitract'), obj.connectome];
                    ADJ_connectome_path = [connectome_folder,filesep,'merged_pathways_ADJ.mat'];
                    obj.ADJ = load(ADJ_connectome_path); % but we need to precompute (ourselves!)
                else
                    [connectome_folder,name,~] = fileparts(cfile);
                    ADJ_connectome_path = [connectome_folder,filesep,name,'_ADJ.mat'];
                end

                try
                    obj.ADJ = load(ADJ_connectome_path);
                catch
                    disp('The adjacency matrix is missing for the chosen connectome')
                    disp('Please, pre-compute it using ea_compute_adj_matrix.m and store the file in the connectome folder')
                    obj.ADJ = false;
                    obj.use_adjacency = false;
                end
            else
                obj.ADJ = false;
            end

            switch obj.connectivity_type
                case 2    % if PAM, then just extracts activation states from fiberActivation.mat
                    pamlist = ea_discfibers_getpams(obj);
                    [fibsvalBin, fibsvalSum, fibsvalMean, fibsvalPeak, fibsval5Peak, fibcell, connFiberInd, totalFibers] = ea_discfibers_calcvals_pam(pamlist, obj, cfile);
                    obj.results.(ea_conn2connid(obj.connectome)).('PAM_Ttest').fibsval = fibsvalBin;
                    obj.results.(ea_conn2connid(obj.connectome)).connFiberInd_PAM = connFiberInd;
                    obj.results.(ea_conn2connid(obj.connectome)).totalFibers = totalFibers; % total number of fibers in the connectome to work with global indices
                otherwise     % check fiber recruitment via intersection with VTA
                    if isfield(obj.M,'pseudoM')
                        vatlist = obj.M.ROI.list;
                    else
                        vatlist = ea_discfibers_getvats(obj);
                    end
                    ea_discfibers_roi_collect(obj); % integrate ROI into .fibfilt file

                    [fibsvalBin, fibsvalSum, fibsvalMean, fibsvalPeak, fibsval5Peak, fibcell,  connFiberInd, totalFibers] = ea_discfibers_calcvals(vatlist, cfile, obj.calcthreshold);
                    obj.results.(ea_conn2connid(obj.connectome)).('VAT_Ttest').fibsval = fibsvalBin;
                    obj.results.(ea_conn2connid(obj.connectome)).connFiberInd_VAT = connFiberInd; % old ff files do not have these data and will fail when using pathway atlases
                    obj.results.(ea_conn2connid(obj.connectome)).totalFibers = totalFibers; % total number of fibers in the connectome to work with global indices
            end

            obj.results.(ea_conn2connid(obj.connectome)).('spearman_sum').fibsval = fibsvalSum;
            obj.results.(ea_conn2connid(obj.connectome)).('spearman_mean').fibsval = fibsvalMean;
            obj.results.(ea_conn2connid(obj.connectome)).('spearman_peak').fibsval = fibsvalPeak;
            obj.results.(ea_conn2connid(obj.connectome)).('spearman_5peak').fibsval = fibsval5Peak;
            obj.results.(ea_conn2connid(obj.connectome)).('plainconn').fibsval = fibsvalBin;
            obj.results.(ea_conn2connid(obj.connectome)).fibcell = fibcell;
        end


        function calculate_cleartune(obj,Efields)
            if isequal(obj.cleartuneefields,Efields) % cleartuneresults already calculated with exact same input.
                return
            else
                obj.cleartuneefields=Efields;

                fibcell=obj.results.(ea_conn2connid(obj.connectome)).fibcell;
                [fibsvalBin, fibsvalSum, fibsvalMean, fibsvalPeak, fibsval5Peak, fibcell] = ea_discfibers_calcvals_cleartune(Efields, fibcell, obj.calcthreshold);

                obj.cleartuneresults.(ea_conn2connid(obj.connectome)).('ttests').fibsval = fibsvalBin;
                obj.cleartuneresults.(ea_conn2connid(obj.connectome)).('VAT_Ttest').fibsval = fibsvalBin; % should be used instead of ttests
                obj.cleartuneresults.(ea_conn2connid(obj.connectome)).('spearman_sum').fibsval = fibsvalSum;
                obj.cleartuneresults.(ea_conn2connid(obj.connectome)).('spearman_mean').fibsval = fibsvalMean;
                obj.cleartuneresults.(ea_conn2connid(obj.connectome)).('spearman_peak').fibsval = fibsvalPeak;
                obj.cleartuneresults.(ea_conn2connid(obj.connectome)).('spearman_5peak').fibsval = fibsval5Peak;
                obj.cleartuneresults.(ea_conn2connid(obj.connectome)).('plainconn').fibsval = fibsvalBin;
                obj.cleartuneresults.(ea_conn2connid(obj.connectome)).fibcell = fibcell;
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
            D = load(obj.leadgroup);
            obj.M = D.M;
            obj.allpatients=obj.M.patient.list;
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
                [I,Ihat, ~, val_struct] = crossval(obj,cvp,[],0,silent);
            else
                % plot some statistics over shuffles
                r_over_iter = zeros(iter,1);
                p_over_iter = zeros(iter,1);
                for i=1:iter
                    cvp = cvpartition(length(obj.patientselection), 'KFold', obj.kfold);
                    if ~silent
                        fprintf("Iterating fold set: %d",i)
                    end
                    [I_iter{i}, Ihat_iter{i}] = crossval(obj, cvp, [], 1);
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
                [I, Ihat, ~, val_struct] = crossval(obj, cvp, [], [], silent);
            else
                [I, Ihat, ~, val_struct] = crossval(obj, cvp, Iperm, [], silent);
            end
        end

        function [Improvement, Ihat, actualimprovs, val_struct] = crossval(obj, cvp, Iperm, shuffle, silent)
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

            if obj.useExternalModel == true && ~strcmp(obj.ExternalModelFile, 'None')
                S = load(obj.ExternalModelFile);
                if ~strcmp(ea_method2methodid(obj),S.fibsvalType)
                    waitfor(msgbox('Change the Model Setup! See terminal'));
                    disp('The loaded model uses: ')
                    disp(S.fibsvalType)
                end

                fibsval = full(obj.results.(ea_conn2connid(obj.connectome)).(S.fibsvalType).fibsval);
            else
                fibsval = full(obj.results.(ea_conn2connid(obj.connectome)).(ea_method2methodid(obj)).fibsval);
            end

            %fibsval = full(obj.results.(ea_conn2connid(obj.connectome)).(ea_method2methodid(obj)).fibsval);

            % for nested LOO, store some statistics
            if obj.nestedLOO
                Abs_pred_error = zeros(cvp.NumTestSets, 1);
                Predicted_scores = zeros(length(patientsel), 1);
                Slope = zeros(cvp.NumTestSets, 1);
                Intercept = zeros(cvp.NumTestSets, 1);

            end

            % no smoothing if the adjacency matrix was not defined
            if ~isstruct(obj.ADJ)
                obj.adj_scaler = 0.0;
            end
            obj.adj_scaler = 0.0;
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
                            [Ihat_inner, ~, ~,actualimprovs] = ea_compute_fibscore_model(c, obj.adj_scaler, obj, fibsval, Ihat_inner, Ihat_train_global_inner, patientsel, training_inner, test_inner);
                        else
                            [Ihat_inner, ~, ~,actualimprovs] = ea_compute_fibscore_model(c, obj.adj_scaler, obj, fibsval, Ihat_inner, Ihat_train_global_inner, patientsel, training_inner, test_inner,Iperm);
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
                    [Ihat, Ihat_train_global, val_struct{c}, actualimprovs] = ea_compute_fibscore_model(c, obj.adj_scaler, obj, fibsval, Ihat, Ihat_train_global, patientsel, training, test);
                else
                    [Ihat, Ihat_train_global, val_struct{c}, actualimprovs] = ea_compute_fibscore_model(c, obj.adj_scaler, obj, fibsval, Ihat, Ihat_train_global, patientsel, training, test, Iperm);
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

            % check if binary variable
            if all(ismember(Improvement(:,1), [0,1])) && size(val_struct{c}.vals,1) == 1
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
                %Ihat_voters_prediction=nan(size(Ihat_voters));
                for c=1:cvp.NumTestSets

                    if isobject(cvp)
                        training = cvp.training(c);
                        test = cvp.test(c);
                    elseif isstruct(cvp)
                        training = cvp.training{c};
                        test = cvp.test{c};
                    end
                    for voter=1:size(val_struct{c}.vals,1)
                        switch obj.multitractmode
                            case 'Split & Color By Subscore'
                                useI=obj.subscore.vars{voter}(patientsel);
                            case 'Split & Color By PCA'
                                useI=obj.subscore.pcavars{voter}(patientsel);
                            otherwise
                                useI=obj.responsevar(patientsel);
                        end

                        if size(useI,2)>1
                            ea_error('This has not been implemented for hemiscores.');
                        end

                        % these predictors are defined within the same ff model
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
                        if ~isempty(covariates)
                            mdl=fitglm([predictor_training(c,training,voter)',covariates(training,:)],useI(training),lower(obj.predictionmodel));
                        else
                            mdl=fitglm([predictor_training(c,training,voter)],useI(training),lower(obj.predictionmodel));
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
                    if obj.CleartuneOptim
                        Ihat = reshape(Ihat,2,length(obj.subscore.vars))';
                    end
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

        function save(obj)
            tractset=obj;
            if isempty(tractset.analysispath)
                pth = fileparts(tractset.leadgroup);
                tractset.analysispath=[pth,filesep,'fiberfiltering',filesep,obj.ID,'.fibfilt'];
                ea_mkdir([pth,filesep,'fiberfiltering']);
            end
            rf=obj.resultfig; % need to stash fig handle for saving.
            rd=obj.drawobject; % need to stash handle of drawing before saving.
            try % could be figure is already closed.
                setappdata(rf,['dt_',tractset.ID],rd); % store handle of tract to figure.
            end

            %we do not need to store plainconn separately since it is a
            %duplicate of PAM or VAT connectivity
            if ~isempty(obj.results)
                if  isfield(obj.results.(ea_conn2connid(obj.connectome)),'plainconn')
                    connectomeName = ea_conn2connid(obj.connectome);
                    obj.results.(connectomeName) = rmfield(obj.results.(connectomeName),'plainconn');
                end
            end

            tractset.useExternalModel = false;
            tractset.ExternalModelFile = 'None'; % do not store imported models
            tractset.resultfig=[]; % rm figure handle before saving.
            tractset.drawobject=[]; % rm drawobject.
            save(tractset.analysispath,'tractset','-v7.3');
            obj.resultfig=rf;
            obj.drawobject=rd;
        end

        function draw(obj,vals,fibcell,usedidx) %for cv live visualize
            %function draw(obj,vals,fibcell)

            % re-define plainconn (since we do not store it)
            try
                if obj.connectivity_type == 2
                    obj.results.(ea_conn2connid(obj.connectome)).('plainconn').fibsval = obj.results.(ea_conn2connid(obj.connectome)).('PAM_Ttest').fibsval;
                else
                    obj.results.(ea_conn2connid(obj.connectome)).('plainconn').fibsval = obj.results.(ea_conn2connid(obj.connectome)).('VAT_Ttest').fibsval;
                end
            catch
                ea_warndlg("Connectivity indices were not stored. Please recalculate or stay with the same model (VAT or PAM)");
                disp("=================== WARNING ========================")
                disp("Connectivity indices connFiberInd were not stored")
                disp("Recalculate or stay with the same model (VAT or PAM)")
                disp("====================================================")
            end

            %disp('Connectivity switch')
            %disp(obj.switch_connectivity)

            % update the connectivity if switched between PAM and VAT
            if obj.switch_connectivity == 1
                if obj.multi_pathways == 1
                    [filepath,name,ext] = fileparts(obj.leadgroup);
                    cfile = [filepath,filesep,'merged_pathways.mat'];
                else
                    cfile = [ea_getconnectomebase('dMRI'), obj.connectome, filesep, 'data.mat'];
                end
                load(cfile, 'fibers', 'idx');
                %disp('Conn. Type:')
                %disp(ea_method2methodid(obj))
                obj.results.(ea_conn2connid(obj.connectome)).totalFibers = length(idx);

                try
                    for side = 1:2
                        if obj.connectivity_type == 2
                            connFiber = fibers(ismember(fibers(:,4), obj.results.(ea_conn2connid(obj.connectome)).connFiberInd_PAM{side}), 1:3);
                            obj.results.(ea_conn2connid(obj.connectome)).fibcell{side} = mat2cell(connFiber, idx(obj.results.(ea_conn2connid(obj.connectome)).connFiberInd_PAM{side}));
                        else
                            connFiber = fibers(ismember(fibers(:,4), obj.results.(ea_conn2connid(obj.connectome)).connFiberInd_VAT{side}), 1:3);
                            obj.results.(ea_conn2connid(obj.connectome)).fibcell{side} = mat2cell(connFiber, idx(obj.results.(ea_conn2connid(obj.connectome)).connFiberInd_VAT{side}));
                        end
                    end
                catch
                    ea_warndlg("Connectivity indices were not stored. Please recalculate or stay with the same model (VAT or PAM)");
                    disp("=================== WARNING ========================")
                    disp("Connectivity indices connFiberInd were not stored")
                    disp("Recalculate or stay with the same model (VAT or PAM)")
                    disp("====================================================")
                end
            else
                % legacy support
                if ~isfield(obj.results.(ea_conn2connid(obj.connectome)),'totalFibers')
                    if obj.multi_pathways == 1
                        [filepath,~,~] = fileparts(obj.leadgroup);
                        cfile = [filepath,filesep,'merged_pathways.mat'];
                    else
                        cfile = [ea_getconnectomebase('dMRI'), obj.connectome, filesep, 'data.mat'];
                    end
                    load(cfile, 'fibers', 'idx');
                    obj.results.(ea_conn2connid(obj.connectome)).totalFibers = length(idx);
                end
            end

            if obj.useExternalModel == true && ~strcmp(obj.ExternalModelFile, 'None')
                S = load(obj.ExternalModelFile);
                if ~strcmp(S.connectome,ea_conn2connid(obj.connectome))
                    waitfor(msgbox('The chosen fibfilt model was computed for another connectome! See terminal'));
                    disp('Model for connectome: ')
                    disp(S.connectome)
                    return
                end

                if obj.connectivity_type ~= S.conn_type
                    waitfor(msgbox('The connectivity methods of imported and current model are different! See terminal'));
                    disp('Connectivity type of imported model: ')
                    disp(obj.connectivity_type)
                    return
                end

                vals_connected = cell(size(S.vals_all,1),size(S.vals_all,2)); % always iterate both sides
                for voter = 1:size(vals_connected,1)
                    for side=1:size(vals_connected,2)
                        try
                            switch obj.connectivity_type
                                case 2
                                    if size(S.vals_all,2) == 2
                                        vals_connected{voter,side} = S.vals_all{voter,side}(obj.results.(ea_conn2connid(obj.connectome)).connFiberInd_PAM{side});
                                    else
                                        vals_connected{voter,side} = S.vals_all{voter,1}(obj.results.(ea_conn2connid(obj.connectome)).connFiberInd_PAM{side});
                                    end
                                otherwise
                                    if size(S.vals_all,2) == 2
                                        vals_connected{voter,side} = S.vals_all{voter,side}(obj.results.(ea_conn2connid(obj.connectome)).connFiberInd_VAT{side});
                                    else
                                        vals_connected{voter,side} = S.vals_all{voter,1}(obj.results.(ea_conn2connid(obj.connectome)).connFiberInd_VAT{side});
                                    end
                            end
                        catch
                            ea_warndlg("Connectivity indices were not stored. Please recalculate.");
                            return
                        end
                    end
                end



                [vals,fibcell,usedidx]=ea_discfibers_loadModel_calcstats(obj, vals_connected);
            elseif ~exist('vals','var')
                [vals,fibcell,usedidx]=ea_discfibers_calcstats(obj);
            end

            ea_discfibers_showroi(obj);

            if ~exist('vals','var')
                [vals,fibcell,usedidx]=ea_discfibers_calcstats(obj);
            end
            obj.fiberdrawn.fibcell = fibcell;
            obj.fiberdrawn.vals = vals;
            obj.fiberdrawn.usedidx = usedidx;

            % if show connected (white) fibers, also calculate those
            if ~isempty(obj.connfiberdrawn)
                for side=1:size(obj.connfiberdrawn.fibcell,2)
                    try
                        delete(obj.conndrawobject{side});
                    end
                end
            end

            if obj.connfibvisible
                obj.runwhite=1;
                [~,cfibcell,~]=ea_discfibers_calcstats(obj);
                obj.runwhite=0;
                obj.connfiberdrawn.fibcell = cfibcell;
                %obj.connfiberdrawn.vals = cvals;
                %obj.connfiberdrawn.usedidx = cusedidx;

                prefs = ea_prefs;
                for side=1:size(obj.connfiberdrawn.fibcell,2)
                    l=length(obj.connfiberdrawn.fibcell{side});
                    % thin out to prefs.fibfilt.maxwhite max
                    if ~isinf(prefs.fibfilt.connfibs.showmax)
                        obj.connfiberdrawn.fibcell{side}=obj.connfiberdrawn.fibcell{side}(1:round(l/prefs.fibfilt.connfibs.showmax):end);
                    end
                    obj.connfiberdrawn.fibcell{side}=ea_discfibers_addjitter(obj.connfiberdrawn.fibcell{side},0.03);


                    obj.conndrawobject{side} = streamtube(obj.connfiberdrawn.fibcell{side}, prefs.fibfilt.connfibs.fiberwidth);
                    nones=repmat({'none'},size(obj.connfiberdrawn.fibcell{side}));
                    [obj.conndrawobject{side}.EdgeColor]=nones{:};

                    % Calulate fiber colors alpha values
                    fibcolor = repmat({prefs.fibfilt.connfibs.color},size(obj.connfiberdrawn.fibcell{side}));
                    fibalpha = repmat({prefs.fibfilt.connfibs.alpha},size(obj.connfiberdrawn.fibcell{side}));

                    % Set fiber colors and alphas
                    [obj.conndrawobject{side}.FaceColor]=fibcolor{:};
                    [obj.conndrawobject{side}.FaceAlpha]=fibalpha{:};
                end
            end

            % print number of significant displayed fibers per pathway (atm only for binary metrics)
            if obj.multi_pathways == 1 && (isequal(ea_method2methodid(obj),'VAT_Ttest') || isequal(ea_method2methodid(obj),'PAM_Ttest') || isequal(ea_method2methodid(obj),'plainconn'))% at the moment, obj.connFiberInd is defined only for OSS-DBS
                %disp("number of drawn fibers per pathway")
                num_per_path = cell(1, 2); % with obj.map_list, rates can be computed
                for side = 1:length(usedidx)
                    num_per_path{side} = zeros(1,length(obj.map_list));
                    if length(usedidx{side})
                        for inx = 1:length(usedidx{side})
                            % check the nearest via the difference, if positive, take one before
                            % I think we can easily add plainconnectivity
                            % here (just try to disable if check)

                            % usedidx will be different for
                            % plainconnectivity, so it should work!
                            if obj.connectivity_type == 2
                                [d, ix] = min(abs(obj.map_list-obj.results.(ea_conn2connid(obj.connectome)).connFiberInd_PAM{side}(usedidx{side}(inx))));
                                if (obj.map_list(ix)-obj.results.(ea_conn2connid(obj.connectome)).connFiberInd_PAM{side}(usedidx{side}(inx))) > 0
                                    ix = ix - 1;
                                end
                            else
                                [d, ix] = min(abs(obj.map_list-obj.results.(ea_conn2connid(obj.connectome)).connFiberInd_VAT{side}(usedidx{side}(inx))));
                                if (obj.map_list(ix)-obj.results.(ea_conn2connid(obj.connectome)).connFiberInd_VAT{side}(usedidx{side}(inx))) > 0
                                    ix = ix - 1;
                                end
                            end
                            num_per_path{side}(ix) = num_per_path{side}(ix)+1;
                        end
                    end
                    %disp(num_per_path{side})  % for now just print number of fibers per pathway
                end

                figure
                t = tiledlayout(1,2,'TileSpacing','compact');
                nonZero_idx = [num_per_path{1}] > 0;
                num_per_path{1} = num_per_path{1}(nonZero_idx);
                if ~isempty(num_per_path{1})
                    % Create pie charts
                    ax1 = nexttile;
                    pie1 = pie(ax1,num_per_path{1});
                    ax1.Colormap = parula(numel(pie1)/2);  % they are all ugl
                    title('Right HS')
                    % Create legend
                    lgd = legend(obj.pathway_list(nonZero_idx));
                    lgd.Layout.Tile = 'west';
                end

                ax2 = nexttile;
                colormap(ax2,winter)
                nonZero_idx = [num_per_path{2}] > 0;
                num_per_path{2} = num_per_path{2}(nonZero_idx);
                if ~isempty(num_per_path{2})
                    pie(ax2,num_per_path{2})
                    title('Left HS')
                    % Create legend
                    lgd2 = legend(obj.pathway_list(nonZero_idx));
                    lgd2.Layout.Tile = 'east';
                end
            end

            allvals{1}=[]; % need to use a loop here - cat doesnt work in all cases with partly empty cells..
            if size(vals,2)==2 % can be a single cell in case of custom code (pseudoM setting).
                allvals{2}=[];
            end
            for v=1:size(vals,1)
                allvals{1}=[allvals{1};vals{v,1}];
                if size(vals,2)==2
                    allvals{2}=[allvals{2};vals{v,2}];
                end
            end
            obj.stats.pos.shown(1)=sum(allvals{1}>0);
            obj.stats.neg.shown(1)=sum(allvals{1}<0);
            if size(vals,2)>1 % bihemispheric usual case
                obj.stats.pos.shown(2)=sum(allvals{2}>0);
                obj.stats.neg.shown(2)=sum(allvals{2}<0);
            end

            set(0,'CurrentFigure',obj.resultfig);

            domultitract=size(vals,1)>1; % if color by groups is set will be positive.
            if ~isfield(obj.M,'groups')
                obj.M.groups.group=1;
                obj.M.groups.color=ea_color_wes('all');
            end

            switch obj.multitractmode
                case 'Split & Color By Group'
                    linecols=obj.M.groups.color;
                case 'Split & Color By Subscore'
                    linecols = obj.subscore.colors;
                case 'Split & Color By PCA'
                    linecols = obj.subscore.pcacolors;
            end

            if isempty(obj.drawobject) % check if prior object has been stored
                obj.drawobject=getappdata(obj.resultfig,['dt_',obj.ID]); % store handle of tract to figure.
            end

            for tract=1:numel(obj.drawobject)
                % try since could run into error when reopening from scratch.
                try
                    delete(obj.drawobject{tract});
                end
            end

            if strcmp(obj.multitractmode,'Single Tract Analysis') || obj.subscore.special_case
                % reset colorbar
                obj.colorbar=[];
                if ~any([obj.posvisible,obj.negvisible])
                    return
                end
            end

            for group=1:size(vals,1) % vals will have 1x2 in case of bipolar drawing and Nx2 in case of group-based drawings (where only positives are shown).
                % vals will also be >1 for subscore tracts
                % Vertcat all values for colorbar construction
                if domultitract && ~obj.subscore.special_case
                    if ~any([obj.subscore.posvisible(group),obj.subscore.negvisible(group)])
                        continue
                    end
                end
                if domultitract
                    obj.subscore.vis.pos_shown(group,1)=sum(vals{group,1}>0);
                    obj.subscore.vis.neg_shown(group,1)=sum(vals{group,1}<0);
                    if size(vals{group},2)==1 % unihemispheric case, should then be pseudoM case
                        obj.subscore.vis.pos_shown(group,2) = 0;
                        obj.subscore.vis.neg_shown(group,2) = 0;
                    elseif (size(vals{group,2},1))>1 % bihemispheric usual case
                        obj.subscore.vis.pos_shown(group,2)=sum(vals{group,2}>0);
                        obj.subscore.vis.neg_shown(group,2)=sum(vals{group,2}<0);
                    elseif length(vals(group,:)) == 2 % in the case that it still exist
                        obj.subscore.vis.pos_shown(group,2) = 0;
                        obj.subscore.vis.neg_shown(group,2) = 0;
                    end
                end
                allvals = full(vertcat(vals{group,:}));

                if isempty(allvals) || all(isnan(allvals))
                    % ea_cprintf('CmdWinWarnings', 'Empty or all-nan value found!\n');
                    continue;
                else
                    allvals(isnan(allvals)) = 0;
                end

                if strcmp(obj.multitractmode,'Split & Color By Subscore') || strcmp(obj.multitractmode,'Split & Color By PCA')
                    if obj.subscore.special_case
                        %basicaly in the mixed fiber case, since the tabs
                        %are communicating with each other, don't set
                        %posfibers & negfibers to be zero unless the vals
                        %across all subscores are pos and neg.
                        %you would also only need to do this once. i.e.,
                        %for group == 1, even if new vals are calculated it
                        %will be taken in account here.
                        if group == 1
                            vals_across_subscores = full(vertcat(vals{1:size(vals,1),:}));
                            if obj.posvisible && all(vals_across_subscores<0)
                                obj.posvisible = 0;
                                fprintf('\n')
                                warning('off', 'backtrace');
                                warning('No positive values found, posvisible is set to 0 now!');
                                warning('on', 'backtrace');
                                fprintf('\n')
                            end

                            if obj.negvisible && all(vals_across_subscores>0)
                                obj.negvisible = 0;
                                fprintf('\n')
                                warning('off', 'backtrace');
                                warning('No negative values found, negvisible is set to 0 now!');
                                warning('on', 'backtrace');
                                fprintf('\n')
                            end
                        end
                    else
                        if obj.subscore.posvisible(group) && all(allvals<0)
                            obj.subscore.posvisible(group) = 0;
                            fprintf('\n')
                            warning('off', 'backtrace');
                            warning('No positive values found, posvisible is set to 0 now!');
                            warning('on', 'backtrace');
                            fprintf('\n')
                        end

                        if obj.subscore.negvisible(group) && all(allvals>0)
                            obj.subscore.negvisible(group) = 0;
                            fprintf('\n')
                            warning('off', 'backtrace');
                            warning('No negative values found, negvisible is set to 0 now!');
                            warning('on', 'backtrace');
                            fprintf('\n')
                        end
                    end

                else
                    if obj.posvisible && all(allvals<0)
                        obj.posvisible = 0;
                        fprintf('\n')
                        warning('off', 'backtrace');
                        warning('No positive values found, posvisible is set to 0 now!');
                        warning('on', 'backtrace');
                        fprintf('\n')
                    end

                    if obj.negvisible && all(allvals>0)
                        obj.negvisible = 0;
                        fprintf('\n')
                        warning('off', 'backtrace');
                        warning('No negative values found, negvisible is set to 0 now!');
                        warning('on', 'backtrace');
                        fprintf('\n')
                    end
                end
                colormap(gray);
                gradientLevel = 1024;
                cmapShiftRatio = 0.5;
                shiftedCmapStart = round(gradientLevel*cmapShiftRatio)+1;
                shiftedCmapEnd = gradientLevel-round(gradientLevel*cmapShiftRatio);
                shiftedCmapLeftEnd = gradientLevel/2-round(gradientLevel/2*cmapShiftRatio);
                shiftedCmapRightStart = round(gradientLevel/2*cmapShiftRatio)+1;

                if domultitract % also means subscores
                    switch obj.multitractmode
                        %logic is different for groups (pos & neg cannot be
                        %shown together), whereas for PCA it is not as
                        %such. Therefore, I am splitting the cases into
                        %two.
                        case 'Split & Color By Group'
                            obj.poscolor = linecols(group,:);
                            obj.negcolor = [0.94,0.97,1.00];
                            if obj.subscore.special_case
                                cmap = ea_colorgradient(gradientLevel, [1,1,1], obj.poscolor);
                                if obj.posvisible && ~obj.negvisible
                                    fibcmap{group} = ea_colorgradient(gradientLevel, cmap(shiftedCmapStart,:), obj.poscolor);
                                    cmapind = round(normalize(allvals,'range',[1,gradientLevel]));
                                    alphaind = ones(size(allvals));
                                    % alphaind = normalize(allvals, 'range');
                                elseif ~obj.posvisible && obj.negvisible
                                    cmap = ea_colorgradient(gradientLevel, obj.negcolor, [1,1,1]);
                                    fibcmap{group} = ea_colorgradient(gradientLevel, obj.negcolor, cmap(shiftedCmapEnd,:));
                                    cmapind = round(normalize(allvals,'range',[1,gradientLevel]));
                                    alphaind = ones(size(allvals));
                                    % alphaind = normalize(-allvals, 'range');
                                else
                                    warndlg(sprintf(['Please choose either "Show Positive Fibers" or "Show Negative Fibers".',...
                                        '\nShow both positive and negative fibers is not supported when "Color by Subscore Variable" is on.']));
                                    return;
                                end
                            else
                                cmap = ea_colorgradient(gradientLevel, [1,1,1], obj.poscolor);
                                if obj.subscore.posvisible(group) && obj.subscore.negvisible(group)
                                    cmap = ea_colorgradient(gradientLevel/2, obj.negcolor, [1,1,1]);
                                    cmapLeft = ea_colorgradient(gradientLevel/2, obj.negcolor, cmap(shiftedCmapLeftEnd,:));
                                    cmap = ea_colorgradient(gradientLevel/2, [1,1,1], obj.poscolor);
                                    cmapRight = ea_colorgradient(gradientLevel/2, cmap(shiftedCmapRightStart,:), obj.poscolor);
                                    fibcmap{group} = [cmapLeft;cmapRight];
                                    cmapind = ones(size(allvals))*gradientLevel/2;
                                    cmapind(allvals<0) = round(normalize(allvals(allvals<0),'range',[1,gradientLevel/2]));
                                    cmapind(allvals>0) = round(normalize(allvals(allvals>0),'range',[gradientLevel/2+1,gradientLevel]));
                                    alphaind = ones(size(allvals));
                                    % alphaind(allvals<0) = normalize(-1./(1+exp(-allvals(allvals<0))), 'range');
                                    % alphaind(allvals>0) = normalize(1./(1+exp(-allvals(allvals>0))), 'range');
                                elseif obj.subscore.posvisible(group) && ~obj.subscore.negvisible(group)
                                    fibcmap{group} = ea_colorgradient(gradientLevel, cmap(shiftedCmapStart,:), obj.poscolor);
                                    cmapind = round(normalize(allvals,'range',[1,gradientLevel]));
                                    alphaind = ones(size(allvals));
                                    % alphaind = normalize(allvals, 'range');
                                elseif ~obj.subscore.posvisible(group) && obj.subscore.negvisible(group)
                                    cmap = ea_colorgradient(gradientLevel, obj.negcolor, [1,1,1]);
                                    fibcmap{group} = ea_colorgradient(gradientLevel, obj.negcolor, cmap(shiftedCmapEnd,:));
                                    cmapind = round(normalize(allvals,'range',[1,gradientLevel]));
                                    alphaind = ones(size(allvals));
                                    % alphaind = normalize(-allvals, 'range');
                                end
                            end
                        otherwise
                            if  strcmp(obj.multitractmode,'Split & Color By Subscore')
                                obj.poscolor = obj.subscore.colors{1,2}(group,:); % positive main color
                                obj.negcolor = obj.subscore.colors{1,1}(group,:); % negative main color
                            elseif strcmp(obj.multitractmode,'Split & Color By PCA')
                                obj.poscolor = obj.subscore.pcacolors(group,:);
                                obj.negcolor = [0.94,0.97,1.00];
                            end
                            if obj.subscore.special_case %operating the mixed fiber tract in the multitract mode. Essentially uses the same logic as you would have used if you were not doing multitract mode, but it incorporates the multitract analysis.
                                %for split by groups options, you cannot have pos &
                                %neg open at the same time.
                                if obj.posvisible && obj.negvisible
                                    cmap = ea_colorgradient(gradientLevel/2, obj.negcolor, [1,1,1]);
                                    cmapLeft = ea_colorgradient(gradientLevel/2, obj.negcolor, cmap(shiftedCmapLeftEnd,:));
                                    cmap = ea_colorgradient(gradientLevel/2, [1,1,1], obj.poscolor);
                                    cmapRight = ea_colorgradient(gradientLevel/2, cmap(shiftedCmapRightStart,:), obj.poscolor);
                                    fibcmap{group} = [cmapLeft;cmapRight];
                                    cmapind = ones(size(allvals))*gradientLevel/2;
                                    cmapind(allvals<0) = round(normalize(allvals(allvals<0),'range',[1,gradientLevel/2]));
                                    cmapind(allvals>0) = round(normalize(allvals(allvals>0),'range',[gradientLevel/2+1,gradientLevel]));
                                    alphaind = ones(size(allvals));
                                    % alphaind(allvals<0) = normalize(-1./(1+exp(-allvals(allvals<0))), 'range');
                                    % alphaind(allvals>0) = normalize(1./(1+exp(-allvals(allvals>0))), 'range');
                                elseif obj.posvisible
                                    cmap = ea_colorgradient(gradientLevel, [1,1,1], obj.poscolor);
                                    fibcmap{group} = ea_colorgradient(gradientLevel, cmap(shiftedCmapStart,:), obj.poscolor);
                                    cmapind = round(normalize(allvals,'range',[1,gradientLevel]));
                                    alphaind = ones(size(allvals));
                                    % alphaind = normalize(1./(1+exp(-allvals)), 'range');
                                elseif obj.negvisible
                                    cmap = ea_colorgradient(gradientLevel, obj.negcolor, [1,1,1]);
                                    fibcmap{group} = ea_colorgradient(gradientLevel, obj.negcolor, cmap(shiftedCmapEnd,:));
                                    cmapind = round(normalize(allvals,'range',[1,gradientLevel]));
                                    alphaind = ones(size(allvals));
                                    % alphaind = normalize(-1./(1+exp(-allvals)), 'range');
                                end
                            else
                                if obj.subscore.posvisible(group) && obj.subscore.negvisible(group)
                                    cmap = ea_colorgradient(gradientLevel/2, obj.negcolor, [1,1,1]);
                                    cmapLeft = ea_colorgradient(gradientLevel/2, obj.negcolor, cmap(shiftedCmapLeftEnd,:));
                                    cmap = ea_colorgradient(gradientLevel/2, [1,1,1], obj.poscolor);
                                    cmapRight = ea_colorgradient(gradientLevel/2, cmap(shiftedCmapRightStart,:), obj.poscolor);
                                    fibcmap{group} = [cmapLeft;cmapRight];
                                    cmapind = ones(size(allvals))*gradientLevel/2;
                                    cmapind(allvals<0) = round(normalize(allvals(allvals<0),'range',[1,gradientLevel/2]));
                                    cmapind(allvals>0) = round(normalize(allvals(allvals>0),'range',[gradientLevel/2+1,gradientLevel]));
                                    alphaind = ones(size(allvals));
                                    % alphaind(allvals<0) = normalize(-1./(1+exp(-allvals(allvals<0))), 'range');
                                    % alphaind(allvals>0) = normalize(1./(1+exp(-allvals(allvals>0))), 'range');
                                elseif obj.subscore.posvisible(group)
                                    cmap = ea_colorgradient(gradientLevel, [1,1,1], obj.poscolor);
                                    fibcmap{group} = ea_colorgradient(gradientLevel, cmap(shiftedCmapStart,:), obj.poscolor);
                                    cmapind = round(normalize(allvals,'range',[1,gradientLevel]));
                                    alphaind = ones(size(allvals));
                                    % alphaind = normalize(1./(1+exp(-allvals)), 'range');
                                elseif obj.subscore.negvisible(group)
                                    cmap = ea_colorgradient(gradientLevel, obj.negcolor, [1,1,1]);
                                    fibcmap{group} = ea_colorgradient(gradientLevel, obj.negcolor, cmap(shiftedCmapEnd,:));
                                    cmapind = round(normalize(allvals,'range',[1,gradientLevel]));
                                    alphaind = ones(size(allvals));
                                    % alphaind = normalize(-1./(1+exp(-allvals)), 'range');
                                end
                            end
                    end
                else
                    obj.poscolor = obj.poscolor; % positive main color
                    obj.negcolor = obj.negcolor; % negative main color

                    if obj.posvisible && obj.negvisible
                        cmap = ea_colorgradient(gradientLevel/2, obj.negcolor, [1,1,1]);
                        cmapLeft = ea_colorgradient(gradientLevel/2, obj.negcolor, cmap(shiftedCmapLeftEnd,:));
                        cmap = ea_colorgradient(gradientLevel/2, [1,1,1], obj.poscolor);
                        cmapRight = ea_colorgradient(gradientLevel/2, cmap(shiftedCmapRightStart,:), obj.poscolor);
                        fibcmap{group} = [cmapLeft;cmapRight];
                        cmapind = ones(size(allvals))*gradientLevel/2;
                        cmapind(allvals<0) = round(normalize(allvals(allvals<0),'range',[1,gradientLevel/2]));
                        cmapind(allvals>0) = round(normalize(allvals(allvals>0),'range',[gradientLevel/2+1,gradientLevel]));
                        alphaind = ones(size(allvals));
                        % alphaind(allvals<0) = normalize(-1./(1+exp(-allvals(allvals<0))), 'range');
                        % alphaind(allvals>0) = normalize(1./(1+exp(-allvals(allvals>0))), 'range');
                    elseif obj.posvisible
                        cmap = ea_colorgradient(gradientLevel, [1,1,1], obj.poscolor);
                        fibcmap{group} = ea_colorgradient(gradientLevel, cmap(shiftedCmapStart,:), obj.poscolor);
                        cmapind = round(normalize(allvals,'range',[1,gradientLevel]));
                        alphaind = ones(size(allvals));
                        % alphaind = normalize(1./(1+exp(-allvals)), 'range');
                    elseif obj.negvisible
                        cmap = ea_colorgradient(gradientLevel, obj.negcolor, [1,1,1]);
                        fibcmap{group} = ea_colorgradient(gradientLevel, obj.negcolor, cmap(shiftedCmapEnd,:));
                        cmapind = round(normalize(allvals,'range',[1,gradientLevel]));
                        alphaind = ones(size(allvals));
                        % alphaind = normalize(-1./(1+exp(-allvals)), 'range');
                    end
                end
                setappdata(obj.resultfig, ['fibcmap',obj.ID], fibcmap);

                if size(vals,2)>1 % standard case
                    cmapind = mat2cell(cmapind, [numel(vals{group,1}), numel(vals{group,2})])';
                    alphaind = mat2cell(alphaind, [numel(vals{group,1}), numel(vals{group,2})])';
                else % potential scripting case, only one side
                    cmapind = mat2cell(cmapind, numel(vals{group,1}))';
                    alphaind = mat2cell(alphaind, numel(vals{group,1}))';
                end
                for side=1:size(vals,2)
                    if domultitract % introduce small jitter for visualization
                        fibcell{group,side}=ea_discfibers_addjitter(fibcell{group,side},0.03);
                    end

                    % Plot fibers if any survived
                    if ~isempty(fibcell{group,side})
                        prefs = ea_prefs;
                        obj.drawobject{group,side} = streamtube(fibcell{group,side}, prefs.d3.fiberwidth);
                        nones=repmat({'none'},size(fibcell{group,side}));
                        [obj.drawobject{group,side}.EdgeColor]=nones{:};

                        % Calulate fiber colors alpha values
                        fibcolor = mat2cell(fibcmap{group}(cmapind{side},:), ones(size(fibcell{group,side})));
                        fibalpha = mat2cell(alphaind{side},ones(size(fibcell{group,side})));

                        % Set fiber colors and alphas
                        [obj.drawobject{group,side}.FaceColor]=fibcolor{:};
                        [obj.drawobject{group,side}.FaceAlpha]=fibalpha{:};
                        obj.drawvals{group,side} = vals{group,side};
                    else
                        obj.drawobject{group,side} = {};
                        obj.drawvals{group,side} = {};
                    end
                end

                % Set colorbar tick positions and labels
                if ~isempty(allvals)
                    if domultitract
                        if obj.subscore.special_case
                            if obj.posvisible && obj.negvisible
                                tick{group} = [1, floor(length(fibcmap{group})/2-40), ceil(length(fibcmap{group})/2+40), length(fibcmap{group})];
                                poscbvals = sort(allvals(allvals>0));
                                negcbvals = sort(allvals(allvals<0));
                                if ~isempty(negcbvals) && ~isempty(poscbvals)
                                    ticklabel{group} = [negcbvals(1), negcbvals(end), poscbvals(1), poscbvals(end)];
                                    ticklabel{group} = arrayfun(@(x) num2str(x,'%.2f'), ticklabel{group}, 'Uni', 0);
                                else
                                    continue
                                end
                            elseif obj.posvisible
                                tick{group} = [1, length(fibcmap{group})];
                                posvals = sort(allvals(allvals>0));
                                ticklabel{group} = [posvals(1), posvals(end)];
                                ticklabel{group} = arrayfun(@(x) num2str(x,'%.2f'), ticklabel{group}, 'Uni', 0);
                            elseif obj.negvisible
                                tick{group} = [1, length(fibcmap{group})];
                                negvals = sort(allvals(allvals<0));
                                ticklabel{group} = [negvals(1), negvals(end)];
                                ticklabel{group} = arrayfun(@(x) num2str(x,'%.2f'), ticklabel{group}, 'Uni', 0);
                            end
                        else
                            if obj.subscore.posvisible(group) && obj.subscore.negvisible(group)
                                tick{group} = [1, floor(length(fibcmap{group})/2-40), ceil(length(fibcmap{group})/2+40), length(fibcmap{group})];
                                poscbvals = sort(allvals(allvals>0));
                                negcbvals = sort(allvals(allvals<0));
                                ticklabel{group} = [negcbvals(1), negcbvals(end), poscbvals(1), poscbvals(end)];
                                ticklabel{group} = arrayfun(@(x) num2str(x,'%.2f'), ticklabel{group}, 'Uni', 0);
                            elseif obj.subscore.posvisible(group)
                                tick{group} = [1, length(fibcmap{group})];
                                posvals = sort(allvals(allvals>0));
                                ticklabel{group} = [posvals(1), posvals(end)];
                                ticklabel{group} = arrayfun(@(x) num2str(x,'%.2f'), ticklabel{group}, 'Uni', 0);
                            elseif obj.subscore.negvisible(group)
                                tick{group} = [1, length(fibcmap{group})];
                                negvals = sort(allvals(allvals<0));
                                ticklabel{group} = [negvals(1), negvals(end)];
                                ticklabel{group} = arrayfun(@(x) num2str(x,'%.2f'), ticklabel{group}, 'Uni', 0);
                            end
                        end
                    else
                        if obj.posvisible && obj.negvisible
                            tick{group} = [1, floor(length(fibcmap{group})/2-40), ceil(length(fibcmap{group})/2+40), length(fibcmap{group})];
                            poscbvals = sort(allvals(allvals>0));
                            negcbvals = sort(allvals(allvals<0));
                            ticklabel{group} = [negcbvals(1), negcbvals(end), poscbvals(1), poscbvals(end)];
                            ticklabel{group} = arrayfun(@(x) num2str(x,'%.2f'), ticklabel{group}, 'Uni', 0);
                        elseif obj.posvisible
                            tick{group} = [1, length(fibcmap{group})];
                            posvals = sort(allvals(allvals>0));
                            ticklabel{group} = [posvals(1), posvals(end)];
                            ticklabel{group} = arrayfun(@(x) num2str(x,'%.2f'), ticklabel{group}, 'Uni', 0);
                        elseif obj.negvisible
                            tick{group} = [1, length(fibcmap{group})];
                            negvals = sort(allvals(allvals<0));
                            ticklabel{group} = [negvals(1), negvals(end)];
                            ticklabel{group} = arrayfun(@(x) num2str(x,'%.2f'), ticklabel{group}, 'Uni', 0);
                        end
                    end
                end
                % store colorbar in object
                if exist('fibcmap','var') % could be no fibers present at all.
                    obj.colorbar.cmap = fibcmap;
                    obj.colorbar.tick = tick;
                    obj.colorbar.ticklabel = ticklabel;
                end

            end
            obj.activate_tractset; % function to highlight tracts activated by a ROI / VTA.
        end

        function activate_tractset(obj)
            if ~isempty(obj.activateby)
                for entry=1:length(obj.activateby)
                    thisentry=obj.activateby{entry};
                    weights={ones(size(obj.cleartuneresults.(ea_conn2connid(obj.connectome)).fibcell{1},1),1),...
                        ones(size(obj.cleartuneresults.(ea_conn2connid(obj.connectome)).fibcell{2},1),1)};
                    if strfind(thisentry,'cleartune')
                        thisentry=strrep(thisentry,'cleartune','');
                        k=strfind(thisentry,'_');
                        ctentry=str2double(thisentry(1:k-1));
                        ctside=str2double(thisentry(k+1:end));
                        weights{ctside}=weights{ctside}+...
                            full(obj.cleartuneresults.(ea_conn2connid(obj.connectome)).(ea_method2methodid(obj)).fibsval{ctside}(:,ctentry));
                    elseif strfind(thisentry,'results')
                        thisentry=strrep(thisentry,'results','');
                        k=strfind(thisentry,'_');
                        ctentry=str2double(thisentry(1:k-1));
                        ctside=str2double(thisentry(k+1:end));
                        weights{ctside}=weights{ctside}+...
                            full(obj.results.(ea_conn2connid(obj.connectome)).(ea_method2methodid(obj)).fibsval{ctside}(:,ctentry));
                    else % manual entry - this could be used to weight a tractset based on a (set of) nifti files.
                        % ignore for now
                    end
                end

                for side=1:size(obj.drawobject,2)
                    if ~(ea_nanmax(weights{side})==1 && ea_nanmin(weights{side})==1)
                        weights{side}=ea_minmax(ea_contrast(weights{side},5))*0.5; % enhance constrast a bit
                    end
                    for entry=1:size(obj.drawobject,1)
                        dweights=weights{side}(obj.fiberdrawn.usedidx{entry,side})';
                        dweights=mat2cell(dweights,1,ones(1,length(dweights)));
                        if ~isempty(dweights)
                            [obj.drawobject{entry,side}.FaceAlpha]=dweights{:};
                        end
                    end
                end
            else
                for side=1:size(obj.drawobject,2)
                    for entry=1:size(obj.drawobject,1)
                        dweights=ones(1,length(obj.drawobject{entry,side}));
                        dweights=mat2cell(dweights,1,ones(1,length(dweights)));
                        try
                            [obj.drawobject{entry,side}.FaceAlpha]=dweights{:};
                        end
                    end
                end
            end
        end
        function copyobj(thisObj,newObj)
            % Construct a new object based on a deep copy of the current
            % object of this class by copying properties over.
            props = properties(thisObj);
            for i = 1:length(props)
                % Use Dynamic Expressions to copy the required property.
                % For more info on usage of Dynamic Expressions, refer to
                % the section "Creating Field Names Dynamically" in:
                % web([docroot '/techdoc/matlab_prog/br04bw6-38.html#br1v5a9-1'])
                newObj.(props{i}) = thisObj.(props{i});
            end
        end

    end

    methods (Static)
        function changeevent(~,event)
            update_trajectory(event.AffectedObject,event.Source.Name);
        end
    end
end


function activatebychange(~,event)
activate_tractset();
end
