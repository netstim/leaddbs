classdef ea_explorerclass < handle
    % Discriminative fiber class to handle visualizations of discriminative fibers in lead dbs resultfig / 3D Matlab figures
    % A. Horn
    properties (SetObservable)
        M % content of lead group project
        resultfig % figure handle to plot results
        ID % name / ID of discriminative fibers object
        resolution = 0.5
        calcthreshold = 50
        modulestate = struct;
        results = struct
        stattests = struct % contains table with information about all statistical tests
        statsettings = struct
        recentmodel = struct
        storedmodels = struct
        thresholding = struct
        multipath = struct        
        %%
        fastrender = 0;
        activateby={}; % entry to use to show fiber activations
        cvlivevisualize = 0; % if set to 1 shows crossvalidation results during processing.
        basepredictionon = 'Mean of Scores';
        fiberdrawn % struct contains fibercell and vals drawn in the resultfig
        drawnstreamlines % actual streamtube handle
        drawnsweetspots
        drawvals % weights of the fibers drawn

        patientselection % selected patients to include. Note that connected fibers are always sampled from all (& mirrored) VTAs of the lead group file
        setlabels={};
        setselections={};
        customselection % selected patients in the custom test list
        allpatients % list of all patients (as from M.patient.list)
        mirrorsides = 1 % flag to mirror VTAs / Efields to contralateral sides using ea_flip_lr_nonlinear()
        responsevar % response variable
        responsevarlabel % label of response variable

        analysispath % where to store results
        leadgroup % redundancy protocol only, path to original lead group project
        connectome % redundancy protocol only, name of underlying connectome
        colorbar % colorbar information

        nestedLOO
        doactualprediction = 0; % set up nested CVs to carry out actual predictions of response variables
        predictionmodel = 'Linear'; % type of glm used to fit fiber values to actual scores
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
        % misc
    end

    properties (Access = private)
        switchedFromSpace=3 % if switching space, this will protocol where from
    end

    methods
        function obj=ea_explorer(analysispath) % class constructor
            if exist('analysispath', 'var') && ~isempty(analysispath)
                obj.analysispath = analysispath;
                [~, ID] = fileparts(obj.analysispath);
                obj.ID = ID;
            end
        end

        function initialize(obj,datapath,resultfig)
            %% set default variables
            obj.modulestate.generation = true;
            obj.modulestate.data = false;
            obj.modulestate.stimstats = false;
            obj.modulestate.threshvis = false;
            obj.modulestate.crosspred = false;
            % Thresholding            
            obj.thresholding.showsignificantonly = 1;
            obj.thresholding.alphalevel = 0.05;
            obj.thresholding.multcompstrategy = 'Uncorrected'; % could be 'Bonferroni' 'FDR'
            obj.thresholding.threshstrategy = 'Percentage Relative to Amount'; % can be 'Relative to Amount' or 'Fixed Amount'
            obj.thresholding.posvisible = 1; % pos tract visible
            obj.thresholding.negvisible = 1; % neg tract visible
            obj.thresholding.showposamount = [100 100]; % two entries for right and left
            obj.thresholding.shownegamount = [100 100]; % two entries for right and left
            obj.thresholding.poscolor = [0.9176,0.2000,0.1373]; % positive main color
            obj.thresholding.negcolor = [0.2824,0.6157,0.9725]; % negative main color
            % List of statistical Tests
            obj.stattests = ea_explorer_statlist;
            % statsettings
            % initial hard threshold to impose on (absolute) nifti files only when calculating the data
            obj.statsettings.doVoxels = 1;
            obj.statsettings.doFibers = 1;
            obj.statsettings.outcometype = 'gradual';
            obj.statsettings.stimulationmodel = 'Electric Field';
            obj.statsettings.efieldmetric = 'Peak'; % if statmetric == ;Correlations / E-fields (Irmen 2020)’, efieldmetric can calculate sum, mean or peak along tracts
            obj.statsettings.efieldthreshold = 200;
            obj.statsettings.connthreshold = 20;
            obj.statsettings.statfamily = 'Correlations'; % the
            obj.statsettings.stattest = 'Spearman';
            obj.statsettings.H0 = 'Average';
            % Multipath
            obj.multipath.multi_pathways = 0; % if structural connectome is devided into pathways (multiple .mat in dMRI_MultiTract)
            obj.multipath.map_list = [];% list that contains global indices of the first fibers in each pathway (relevant when multi_pathways = 1)
            obj.multipath.pathway_list = []; % list that contains names of pathways (relevant when multi_pathways = 1
            %% get data
            datapath = GetFullPath(datapath);
            D = load(datapath, '-mat');
            if isfield(D, 'M') % Lead Group analysis path loaded
                obj.M = D.M;
                obj.leadgroup = datapath;

                testID = obj.M.guid;
                ea_mkdir([fileparts(obj.leadgroup),filesep,'explorer',filesep]);
                id = 1;
                while exist([fileparts(obj.leadgroup),filesep,'explorer',filesep,testID,'.explorer'],'file')
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
                    obj.allpatients = obj.M.patient.list;
                    obj.patientselection = 1:numel(obj.M.patient.list); % always initialize with all patients
                end
                obj.responsevarlabel = obj.M.clinical.labels{1}; % start with first clinical variable
                obj.responsevar = obj.M.clinical.vars{1};
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
            % obj.compat_statmetric; % check and resolve for old statmetric code (which used to be integers)

            addlistener(obj,'activateby','PostSet',@activatebychange);

            % added a check here otherwise errors out for files w/o vatmodels
            if ~isfield(obj.M,'pseudoM')
                if ~isempty(obj.M.vatmodel) && contains(obj.M.vatmodel, 'OSS-DBS (Butenko 2020)')
                    obj.statmetric = 3;
                end
            end
        end

        %% This function gathers Efields and calculates fiber activations
        function calculate(obj)
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
            if obj.multipath.multi_pathways == 1
                [cfile, obj.multipath.map_list, obj.multipath.pathway_list] = ea_discfibers_merge_pathways(obj);
            else
                cfile = [ea_getconnectomebase('dMRI'), obj.connectome, filesep, 'data.mat'];
            end
            % if multi_pathway = 1, use the adjacency defined for
            % merged_pathways.mat

            if isfield(obj.M,'pseudoM')
                vatlist = obj.M.ROI.list;
            else
                vatlist = ea_discfibers_getvats(obj);
            end
            %% Calc Voxels
            [AllX,space] = ea_explorer_calcvoxelvals(vatlist,obj);
            obj.results.efield = AllX;
            obj.results.space = space;
            %% Calc Fibers
            [fibsvalBin, fibsvalSum, fibsvalMean, fibsvalPeak, ~, fibcell,  connFiberInd, totalFibers] = ea_explorer_calcfibervals(vatlist, cfile, obj.calcthreshold);
            obj.results.(ea_conn2connid(obj.connectome)).connFiberInd_VAT = connFiberInd; % old ff files do not have these data and will fail when using pathway atlases
            obj.results.(ea_conn2connid(obj.connectome)).totalFibers = totalFibers; % total number of fibers in the connectome to work with global indices

            obj.results.(ea_conn2connid(obj.connectome)).('sum').fibsval = fibsvalSum;
            obj.results.(ea_conn2connid(obj.connectome)).('mean').fibsval = fibsvalMean;
            obj.results.(ea_conn2connid(obj.connectome)).('peak').fibsval = fibsvalPeak;
            % obj.results.(ea_conn2connid(obj.connectome)).('peak5').fibsval = fibsval5Peak;
            obj.results.(ea_conn2connid(obj.connectome)).('VAT_Ttest').fibsval = fibsvalBin;
            obj.results.(ea_conn2connid(obj.connectome)).fibcell = fibcell;
        end
        %% This function recalculates the main statistical analysis
        function recalculate(obj) %for cv live visualize      
            warning('on','all')
            % re-define plainconn (since we do not store it)
            obj.results.(ea_conn2connid(obj.connectome)).('plainconn').fibsval = obj.results.(ea_conn2connid(obj.connectome)).('VAT_Ttest').fibsval;

            obj.recentmodel = struct;
            obj.recentmodel.resolution = obj.resolution;
            obj.recentmodel.calcthreshold = obj.calcthreshold;
            obj.recentmodel.patientselection = obj.patientselection;
            obj.recentmodel.responsevarlabel = obj.responsevarlabel;
            obj.recentmodel.mirrorsides = obj.mirrorsides;
            obj.recentmodel.statsettings = obj.statsettings;
            obj.recentmodel.thresholding = obj.thresholding;
            if obj.recentmodel.statsettings.doVoxels
                disp('Calculating Voxel Statistics')
                [obj.recentmodel.voxels.vals,obj.recentmodel.voxels.pvals]=ea_explorer_calcstats(obj,'voxels');
            end
            if obj.recentmodel.statsettings.doFibers
                disp('Calculating Fiber Statistics')
                [obj.recentmodel.fibers.vals,obj.recentmodel.fibers.pvals]=ea_explorer_calcstats(obj,'fibers');            
            end
        end
        %% This function draws statistical results like fibers and Sweetspots
        function visualize(obj)
            obj.recentmodel.thresholding = obj.thresholding;
            if obj.recentmodel.statsettings.doVoxels
                disp('Visualizing Voxels')
                ea_explorer_visualizevoxels(obj);
            end
            if obj.recentmodel.statsettings.doFibers
                disp('Visualizing Fibers')
                ea_explorer_visualizefibers(obj);
            end
        end
        %% Whatever these other functions are
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
                    [I_iter{i}, Ihat_iter{i},~,val_struct] = crossval(obj, cvp, [], 1,silent);
                    if ~silent
                        inx_nnan = find(isnan(I_iter{i}) ~= 1);
                        [r_over_iter(i),p_over_iter(i)]=ea_permcorr(I_iter{i}(inx_nnan),Ihat_iter{i}(inx_nnan),'spearman');

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


            if ~exist('Iperm', 'var') || isempty(Iperm)
                Improvement = obj.responsevar(patientsel,:);
            else
                Improvement = Iperm(patientsel,:);
            end


            % Ihat is the estimate of improvements (not scaled to real improvements)

            Ihat = nan(length(patientsel),2);
            Ihat_train_global = nan(cvp.NumTestSets,length(patientsel),2);


            fibsval = full(obj.results.(ea_conn2connid(obj.connectome)).(ea_explorer_method2methodid(obj)).fibsval);


            %fibsval = full(obj.results.(ea_conn2connid(obj.connectome)).(ea_explorer_method2methodid(obj)).fibsval);

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




            Ihat=squeeze(Ihat);
            %Ihat=squeeze(Ihat_voters);

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
                obj.recalculate;
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
                tractset.analysispath=[pth,filesep,'explorer',filesep,obj.ID,'.explorer'];
                ea_mkdir([pth,filesep,'explorer']);
            end
            rf=obj.resultfig; % need to stash fig handle for saving.
            rd=obj.drawnstreamlines; % need to stash handle of drawing before saving.
            try % could be figure is already closed.
                setappdata(rf,['dt_',tractset.ID],rd); % store handle of tract to figure.
            end

            % tractset.useExternalModel = false;
            % tractset.ExternalModelFile = 'None'; % do not store imported models
            tractset.resultfig=[]; % rm figure handle before saving.
            tractset.drawnstreamlines=[]; % rm drawnstreamlines.
            save(tractset.analysispath,'tractset','-v7.3');
            obj.resultfig=rf;
            obj.drawnstreamlines=rd;
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