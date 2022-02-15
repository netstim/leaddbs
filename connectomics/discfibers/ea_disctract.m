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
        showposamount = [25 25] % two entries for right and left
        shownegamount = [25 25] % two entries for right and left
        connthreshold = 20
        efieldthreshold = 200
        statmetric = 1 % stats metric to use, 1 = ttest, 2 = correlations, 3 = OSS DBS pathway activations, 4 = Dice Coeff / VTAs for binary variables, 5 = reverse t-tests & e-fields for binary variables, 6 = show plain connections (no stats)
        multi_pathways = 0 % if structural connectome is devided into pathways (multiple .mat in dMRI_MultiTract)
        map_list % list that contains global indices of the first fibers in each pathway (relevant when multi_pathways = 1) 
        pathway_list % list that contains names of pathways (relevant when multi_pathways = 1
        connFiberInd % list of indices of activated (connected) fibers (relevant when multi_pathways = 1)
        corrtype = 'Spearman' % correlation strategy in case of statmetric == 2.
        efieldmetric = 'Peak' % if statmetric == 2, efieldmetric can calculate sum, mean or peak along tracts
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
        activateby={}; % entry to use to show fiber activations
        cvlivevisualize = 0; % if set to 1 shows crossvalidation results during processing.
        basepredictionon = 'Mean of Scores';
        fiberdrawn % struct contains fibercell and vals drawn in the resultfig
        drawobject % actual streamtube handle
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
        groupcolors = ea_color_wes('all');
        % stats: (how many fibers available and shown etc for GUI)
        modelNormalization='None';
        numBins=15;
        stats
        % additional settings:
        rngseed = 'default';
        Nperm = 1000 % how many permutations in leave-nothing-out permtest strategy
        kfold = 5 % divide into k sets when doing k-fold CV
        Nsets = 5 % divide into N sets when doing Custom (random) set test
        adjustforgroups = 1 % adjust correlations for group effects
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
                    obj.M.patient.list=obj.M.ROI.list; % copies
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
                obj.subscore.spitbysubscore = 0;
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

             addlistener(obj,'activateby','PostSet',...
                @activatebychange);
            %added a check here otherwise errors out for files w/o
            %vatmodels
            if ~isempty(obj.M.vatmodel) && contains(obj.M.vatmodel, 'OSS-DBS (Butenko 2020)')
                obj.statmetric = 3;
            end
            
            % just for now, ask Nanditha to add dMRI_MultiTract files as
            % option in the GUI
            %obj.connectome = 
            
        end


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
            if obj.multi_pathways == 1
                [cfile, obj.map_list, obj.pathway_list] = ea_discfibers_merge_pathways(obj);
            else
                cfile = [ea_getconnectomebase('dMRI'), obj.connectome, filesep, 'data.mat'];
            end
            switch obj.statmetric
                case 3    % if PAM, then just extracts activation states from fiberActivation.mat
                    pamlist = ea_discfibers_getpams(obj);
                    [fibsvalBin, fibsvalSum, fibsvalMean, fibsvalPeak, fibsval5Peak, fibcell, obj.connFiberInd] = ea_discfibers_calcvals_pam(pamlist, obj, cfile);
                otherwise                    
                    % this is not needed for OSS-DBS
                    if isfield(obj.M,'pseudoM')
                        vatlist = obj.M.ROI.list;
                    else
                        vatlist = ea_discfibers_getvats(obj);
                    end
                    [fibsvalBin, fibsvalSum, fibsvalMean, fibsvalPeak, fibsval5Peak, fibcell, obj.connFiberInd] = ea_discfibers_calcvals(vatlist, cfile, obj.calcthreshold);
            end
            
            obj.results.(ea_conn2connid(obj.connectome)).('ttests').fibsval = fibsvalBin;
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

        function [I, Ihat] = loocv(obj)
            rng(obj.rngseed);
            cvp = cvpartition(length(obj.patientselection), 'LeaveOut');
            [I, Ihat] = crossval(obj, cvp);
        end

        function [I, Ihat] = lococv(obj)
            if length(unique(obj.M.patient.group(obj.patientselection))) == 1
                ea_error(sprintf(['Only one cohort in the analysis.\n', ...
                    'Leave-One-Cohort-Out cross-validation not possible.']));
            end
            [I, Ihat] = crossval(obj, obj.M.patient.group(obj.patientselection));
        end

        function [I, Ihat] = kfoldcv(obj)
            rng(obj.rngseed);
            cvp = cvpartition(length(obj.patientselection), 'KFold', obj.kfold);
            [I, Ihat] = crossval(obj, cvp);
        end

        function [I, Ihat] = lno(obj, Iperm)
            rng(obj.rngseed);
            cvp = cvpartition(length(obj.patientselection), 'resubstitution');
            if ~exist('Iperm', 'var')
                [I, Ihat] = crossval(obj, cvp);
            else
                [I, Ihat] = crossval(obj, cvp, Iperm);
            end
        end

        function [Improvement, Ihat, actualimprovs] = crossval(obj, cvp, Iperm)
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
                        Improvement = obj.subscore.vars;
                    else
                        ea_error('Permutation based test not yet coded in for PCA mode.');
                    end
                otherwise
                    if ~exist('Iperm', 'var') || isempty(Iperm)
                        Improvement = obj.responsevar(patientsel,:);
                    else
                        Improvement = Iperm(patientsel,:);
                    end
            end
            % Ihat is the estimate of improvements (not scaled to real improvements)

            Ihat = nan(length(patientsel),2);
            Ihattrain = Ihat;

            fibsval = full(obj.results.(ea_conn2connid(obj.connectome)).(ea_method2methodid(obj)).fibsval);

            for c=1:cvp.NumTestSets
                if cvp.NumTestSets ~= 1
                    fprintf(['\nIterating set: %0',num2str(numel(num2str(cvp.NumTestSets))),'d/%d\n'], c, cvp.NumTestSets);
                end

                if isobject(cvp)
                    training = cvp.training(c);
                    test = cvp.test(c);
                elseif isstruct(cvp)
                    training = cvp.training{c};
                    test = cvp.test{c};
                end

                if ~exist('Iperm', 'var')
                    if obj.cvlivevisualize
                        [vals,fibcell,usedidx] = ea_discfibers_calcstats(obj, patientsel(training));
                        obj.draw(vals,fibcell,usedidx)
                        %obj.draw(vals,fibcell);
                        drawnow;
                    else
                        [vals,~,usedidx] = ea_discfibers_calcstats(obj, patientsel(training));
                    end
                else
                    if obj.cvlivevisualize
                        [vals,fibcell,usedidx] = ea_discfibers_calcstats(obj, patientsel(training), Iperm);
                        obj.draw(vals,fibcell,usedidx)
                        %obj.draw(vals,fibcell);
                        drawnow;
                    else
                        [vals,~,usedidx] = ea_discfibers_calcstats(obj, patientsel(training), Iperm);
                    end
                end

                switch obj.modelNormalization
                    case 'z-score'
                        for s=1:length(vals)
                            vals{s}=ea_nanzscore(vals{s});
                        end
                    case 'van Albada 2007'
                        for s=1:length(vals)
                            vals{s}=ea_normal(vals{s});
                        end
                end
                Ihat_voters=[];
                for voter=1:size(vals,1)
                    for side=1:size(vals,2)
                        if ~isempty(vals{voter,side})
                            switch obj.statmetric % also differentiate between methods in the prediction part.
                                case {1,3,4} % ttests / OSS-DBS / reverse t-tests
                                    switch lower(obj.basepredictionon)
                                        case 'mean of scores'
                                            Ihat(test,side) = ea_nanmean(vals{voter,side}.*fibsval{1,side}(usedidx{voter,side},patientsel(test)),1);
                                            if nargout>2
                                                Ihattrain(training,side) = ea_nanmean(vals{voter,side}.*fibsval{1,side}(usedidx{voter,side},patientsel(training)),1);
                                            end
                                        case 'sum of scores'
                                            Ihat(test,side) = ea_nansum(vals{voter,side}.*fibsval{1,side}(usedidx{voter,side},patientsel(test)),1);
                                            if nargout>2
                                                Ihattrain(training,side) = ea_nansum(vals{voter,side}.*fibsval{1,side}(usedidx{voter,side},patientsel(training)),1);
                                            end
                                        case 'peak of scores'
                                            Ihat(test,side) = ea_nanmax(vals{voter,side}.*fibsval{1,side}(usedidx{voter,side},patientsel(test)),1);
                                            if nargout>2
                                                Ihattrain(training,side) = ea_nanmax(vals{voter,side}.*fibsval{1,side}(usedidx{voter,side},patientsel(training)),1);
                                            end
                                        case 'peak 5% of scores'
                                            ihatvals=vals{1,side}.*fibsval{1,side}(usedidx{voter,side},patientsel(test));
                                            ihatvals=sort(ihatvals);
                                            Ihat(test,side) = ea_nansum(ihatvals(1:ceil(size(ihatvals,1).*0.05),:),1);
                                             if nargout>2
                                                ihatvals=vals{1,side}.*fibsval{1,side}(usedidx{voter,side},patientsel(training));
                                                ihatvals=sort(ihatvals);
                                                Ihattrain(training,side) = ea_nansum(ihatvals(1:ceil(size(ihatvals,1).*0.05),:),1);
                                            end

                                    end
                                case {2,5} % efields
                                    switch lower(obj.basepredictionon)
                                        case 'profile of scores: spearman'
                                            Ihat(test,side) = atanh(corr(vals{voter,side},fibsval{1,side}(usedidx{voter,side},patientsel(test)),'rows','pairwise','type','spearman'));
                                            if any(isnan(Ihat(test,side)))
                                                Ihat(isnan(Ihat(test,side)),side)=0;
                                                warning('Profiles of scores could not be evaluated for some patients. Displaying these points as zero entries. Lower threshold or carefully check results.');
                                            end
                                            if nargout>2
                                                Ihattrain(training,side) = atanh(corr(vals{voter,side},fibsval{1,side}(usedidx{voter,side},patientsel(training)),'rows','pairwise','type','spearman'));
                                            end
                                        case 'profile of scores: pearson'
                                            Ihat(test,side) = atanh(corr(vals{voter,side},fibsval{1,side}(usedidx{voter,side},patientsel(test)),'rows','pairwise','type','pearson'));
                                            if any(isnan(Ihat(test,side)))
                                                Ihat(isnan(Ihat(test,side)),side)=0;
                                                warning('Profiles of scores could not be evaluated for some patients. Displaying these points as zero entries. Lower threshold or carefully check results.');
                                            end
                                            if nargout>2
                                                Ihattrain(training,side) = atanh(corr(vals{voter,side},fibsval{1,side}(usedidx{voter,side},patientsel(training)),'rows','pairwise','type','pearson'));
                                            end
                                        case 'profile of scores: bend'
                                            Ihat(test,side) = atanh(ea_bendcorr(vals{voter,side},fibsval{1,side}(usedidx{voter,side},patientsel(test))));
                                            if any(isnan(Ihat(test,side)))
                                                Ihat(isnan(Ihat(test,side)),side)=0;
                                                warning('Profiles of scores could not be evaluated for some patients. Displaying these points as zero entries. Lower threshold or carefully check results.');
                                            end
                                            if nargout>2
                                                Ihattrain(training,side) = atanh(ea_bendcorr(vals{voter,side},fibsval{1,side}(usedidx{voter,side},patientsel(training))));
                                            end
                                        case 'mean of scores'
                                            if ~isempty(vals{voter,side})
                                                Ihat(test,side) = ea_nanmean(vals{voter,side}.*fibsval{1,side}(usedidx{voter,side},patientsel(test)),1);
                                            end
                                            if nargout>2
                                                Ihattrain(training,side) = ea_nanmean(vals{voter,side}.*fibsval{1,side}(usedidx{voter,side},patientsel(training)),1);
                                            end
                                        case 'sum of scores'
                                            if ~isempty(vals{voter,side})
                                                Ihat(test,side) = ea_nansum(vals{voter,side}.*fibsval{1,side}(usedidx{voter,side},patientsel(test)),1);
                                            end
                                            if nargout>2
                                                Ihattrain(training,side) = ea_nansum(vals{voter,side}.*fibsval{1,side}(usedidx{voter,side},patientsel(training)),1);
                                            end
                                        case 'peak of scores'
                                            if ~isempty(vals{voter,side})
                                                Ihat(test,side) = ea_nanmax(vals{voter,side}.*fibsval{1,side}(usedidx{voter,side},patientsel(test)),1);
                                            end
                                            if nargout>2
                                                Ihattrain(training,side) = ea_nanmax(vals{voter,side}.*fibsval{1,side}(usedidx{voter,side},patientsel(training)),1);
                                            end
                                        case 'peak 5% of scores'
                                            if ~isempty(vals{voter,side})
                                                ihatvals=vals{voter,side}.*fibsval{1,side}(usedidx{voter,side},patientsel(test));
                                            end
                                            ihatvals=sort(ihatvals);
                                            Ihat(test,side) = ea_nansum(ihatvals(1:ceil(size(ihatvals,1).*0.05),:),1);
                                            if nargout>2
                                                ihatvals=vals{voter,side}.*fibsval{1,side}(usedidx{voter,side},patientsel(training));
                                                ihatvals=sort(ihatvals);
                                                Ihattrain(training,side) = ea_nansum(ihatvals(1:ceil(size(ihatvals,1).*0.05),:),1);
                                            end
                                    end

                                case 6 % Plain Connection
                                    switch lower(obj.basepredictionon)
                                        case 'mean of scores'
                                            Ihat(test,side) = ea_nanmean(vals{voter,side}.*fibsval{1,side}(usedidx{voter,side},patientsel(test)),1);
                                            if nargout>2
                                                Ihattrain(training,side) = ea_nanmean(vals{voter,side}.*fibsval{1,side}(usedidx{voter,side},patientsel(training)),1);
                                            end
                                        case 'sum of scores'
                                            Ihat(test,side) = ea_nansum(vals{voter,side}.*fibsval{1,side}(usedidx{voter,side},patientsel(test)),1);
                                            if nargout>2
                                                Ihattrain(training,side) = ea_nansum(vals{voter,side}.*fibsval{1,side}(usedidx{voter,side},patientsel(training)),1);
                                            end
                                        case 'peak of scores'
                                            Ihat(test,side) = ea_nanmax(vals{voter,side}.*fibsval{1,side}(usedidx{voter,side},patientsel(test)),1);
                                            if nargout>2
                                                Ihattrain(training,side) = ea_nanmax(vals{voter,side}.*fibsval{1,side}(usedidx{voter,side},patientsel(training)),1);
                                            end
                                        case 'peak 5% of scores'
                                            ihatvals=vals{voter,side}.*fibsval{1,side}(usedidx{voter,side},patientsel(test));
                                            ihatvals=sort(ihatvals);
                                            Ihat(test,side) = ea_nansum(ihatvals(1:ceil(size(ihatvals,1).*0.05),:),1);
                                            if nargout>2
                                                ihatvals=vals{voter,side}.*fibsval{1,side}(usedidx{voter,side},patientsel(training));
                                                ihatvals=sort(ihatvals);
                                                Ihattrain(training,side) = ea_nansum(ihatvals(1:ceil(size(ihatvals,1).*0.05),:),1);
                                            end
                                    end
                            end
                        end
                    end

                    if nargout>2 % send out improvements of subscores

                        switch obj.multitractmode
                            case 'Split & Color By Subscore'
                                useI=obj.subscore.vars{voter};
                            case 'Split & Color By PCA'
                                useI=obj.subscore.pcavars{voter};
                            otherwise
                                useI=obj.responsevar;
                        end
                        for side=1:2
                            mdl=fitglm(Ihattrain(training,side),useI(training),lower(obj.predictionmodel));
                            actualimprovs{voter,side}=predict(mdl,Ihat(test,side));
                        end

                    end

                    Ihat_voters=cat(3,Ihat_voters,Ihat);
                end
            end

            if obj.doactualprediction % repeat loops partly to fit to actual response variables:
                Ihat_voters_prediction=nan(size(Ihat_voters));
                for c=1:cvp.NumTestSets

                    if isobject(cvp)
                        training = cvp.training(c);
                        test = cvp.test(c);
                    elseif isstruct(cvp)
                        training = cvp.training{c};
                        test = cvp.test{c};
                    end
                    for voter=1:size(vals,1)
                        switch obj.multitractmode
                            case 'Split & Color By Subscore'
                                useI=obj.subscore.vars{voter};
                            case 'Split & Color By PCA'
                                useI=obj.subscore.pcavars{voter};
                            otherwise
                                useI=obj.responsevar;
                        end

                        if size(useI,2)>1
                            ea_error('This has not been implemented for hemiscores.');
                        end
                        predictor=squeeze(ea_nanmean(Ihat_voters,2));

                        mdl=fitglm(predictor(training,voter),useI(training),lower(obj.predictionmodel));

                        if size(useI,2)==1 % global scores
                            Ihat_voters_prediction(test,:,voter)=repmat(predict(mdl,predictor(test,voter)),1,2); % fill both sides equally
                        end
                    end
                end
                Ihat_voters=Ihat_voters_prediction; % replace with actual response variables.
            end

            switch obj.multitractmode
                case 'Split & Color By Subscore'
                    % here we map back to the single response variable using a
                    % weightmatrix
                    if isempty(obj.customselection)
                        selected_pts = obj.patientselection;
                    else
                        selected_pts = obj.customselection;
                    end
                    weightmatrix=zeros(size(Ihat_voters));
                    for voter=1:size(Ihat_voters,3)
                        if ~isnan(obj.subscore.weights(voter)) % same weight for all subjects in that voter (slider was used)
                            weightmatrix(:,:,voter)=obj.subscore.weights(voter);
                        else % if the weight value is nan, this means we will need to derive a weight from the variable of choice
                            weightmatrix(:,:,voter)=repmat(obj.subscore.weightvars{voter}(selected_pts),1,size(weightmatrix,2)/size(obj.subscore.weightvars{voter}(selected_pts),2));
                        end
                    end
                    for xx=1:size(Ihat_voters,1) % make sure voter weights sum up to 1
                        for yy=1:size(Ihat_voters,2)
                            weightmatrix(xx,yy,:)=weightmatrix(xx,yy,:)./ea_nansum(weightmatrix(xx,yy,:));
                        end
                    end

                    Ihat=ea_nansum(Ihat_voters.*weightmatrix,3);
                case 'Split & Color By PCA'

                    Ihat_voters=squeeze(ea_nanmean(Ihat_voters,2)); % need to assume global scores here for now.

                    % map back to PCA:
                    subvars=ea_nanzscore(cell2mat(obj.subscore.vars'));
                    [coeff,score,latent,tsquared,explained,mu]=pca(subvars,'rows','pairwise');

                    Ihatout = Ihat_voters*coeff(:,1:obj.numpcs)' + repmat(mu,size(score,1),1);

                    Ihat=cell(1); % export in cell format as the Improvement itself.
                    for subsc=1:size(Ihatout,2)
                        Ihat{subsc}=Ihatout(:,subsc);
                    end

                otherwise
                    Ihat=squeeze(Ihat_voters);
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


        function [Iperm, Ihat, R0, R1, pperm, Rp95] = lnopb(obj, corrType)
            if ~exist('corrType', 'var')
                corrType = 'Spearman';
            end

            numPerm = obj.Nperm;

            Iperm = ea_shuffle(obj.responsevar, numPerm, obj.patientselection, obj.rngseed)';
            Iperm = [obj.responsevar, Iperm];
            Ihat = cell(numPerm+1, 1);

            R = zeros(numPerm+1, 1);

            for perm=1:numPerm+1
                if perm==1
                    fprintf('Calculating without permutation\n\n');
                    [~, Ihat{perm}] = lno(obj);
                else
                    fprintf('Calculating permutation: %d/%d\n\n', perm-1, numPerm);
                    [~, Ihat{perm}] = lno(obj, Iperm(:, perm));
                end

                R(perm) = corr(Iperm(obj.patientselection,perm),Ihat{perm},'type',corrType,'rows','pairwise');
            end

            % generate null distribution
            R1 = R(1);
            R0 = sort((R(2:end)),'descend');

            Rp95 = R0(round(0.05*numPerm));
            higherR0=R0>(R1);
            pperm=sum(higherR0)/numel(higherR0);

            % old method below, using the one above to be consistent with
            % ea_plothistperm.

            %v = ea_searchclosest(R0, R1);
            %pperm = v/numPerm;
            disp(['Permuted p = ',sprintf('%0.2f',pperm),'.']);

            % Return only selected I
            Iperm = Iperm(obj.patientselection,:);
        end

        function save(obj)
            tractset=obj;
            pth = fileparts(tractset.leadgroup);
            tractset.analysispath=[pth,filesep,'fiberfiltering',filesep,obj.ID,'.fibfilt'];
            ea_mkdir([pth,filesep,'fiberfiltering']);
            rf=obj.resultfig; % need to stash fig handle for saving.
            rd=obj.drawobject; % need to stash handle of drawing before saving.
            try % could be figure is already closed.
                setappdata(rf,['dt_',tractset.ID],rd); % store handle of tract to figure.
            end
            tractset.resultfig=[]; % rm figure handle before saving.
            tractset.drawobject=[]; % rm drawobject.
            save(tractset.analysispath,'tractset','-v7.3');
            obj.resultfig=rf;
            obj.drawobject=rd;
        end
        
        function draw(obj,vals,fibcell,usedidx) %for cv live visualize
        %function draw(obj,vals,fibcell)
            if ~exist('vals','var')
                [vals,fibcell,usedidx]=ea_discfibers_calcstats(obj);
            end
            obj.fiberdrawn.fibcell = fibcell;
            obj.fiberdrawn.vals = vals;
            obj.fiberdrawn.usedidx = usedidx;
            
            % print number of significant displayed fibers per pathway
            if obj.multi_pathways == 1 % at the moment, obj.connFiberInd is defined only for OSS-DBS
                %disp("number of drawn fibers per pathway")
                num_per_path = cell(1, 2); % with obj.map_list, rates can be computed                
                for side = 1:2
                    num_per_path{side} = zeros(1,length(obj.map_list));
                    if length(usedidx{side}) 
                        for inx = 1:length(usedidx{side})
                            % check the nearest via the difference, if positive, take one before
                            [d, ix] = min(abs(obj.map_list-obj.connFiberInd{side}(usedidx{side}(inx)))); 
                            if (obj.map_list(ix)-obj.connFiberInd{side}(usedidx{side}(inx))) > 0
                               ix = ix - 1; 
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
                    if (size(vals{group,2},1))>1 % bihemispheric usual case
                        obj.subscore.vis.pos_shown(group,2)=sum(vals{group,2}>0);
                        obj.subscore.vis.neg_shown(group,2)=sum(vals{group,2}<0);
                    elseif length(vals(group,:)) == 2 %in the case that it still exist
                        obj.subscore.vis.pos_shown(group,2) = 0;
                        obj.subscore.vis.neg_shown(group,2) = 0;
                    end
                end
                allvals = full(vertcat(vals{group,:}));

                if isempty(allvals)
                    continue;
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
                            obj.poscolor = obj.groupcolors(group,:);
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
                                if obj.subscore.posvisible(group) && ~obj.subscore.negvisible(group)
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
                                else
                                    warndlg(sprintf(['Please choose either "Show Positive Fibers" or "Show Negative Fibers".',...
                                        '\nShow both positive and negative fibers is not supported when "Color by Subscore Variable" is on.']));
                                    return;
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
                    obj.poscolor = [0.9176,0.2000,0.1373]; % positive main color
                    obj.negcolor = [0.2824,0.6157,0.9725]; % negative main color

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
                        obj.drawobject{group,side}=streamtube(fibcell{group,side},0.2);
                        nones=repmat({'none'},size(fibcell{group,side}));
                        [obj.drawobject{group,side}.EdgeColor]=nones{:};

                        % Calulate fiber colors alpha values
                        fibcolor = mat2cell(fibcmap{group}(cmapind{side},:), ones(size(fibcell{group,side})));
                        fibalpha = mat2cell(alphaind{side},ones(size(fibcell{group,side})));

                        % Set fiber colors and alphas
                        [obj.drawobject{group,side}.FaceColor]=fibcolor{:};
                        [obj.drawobject{group,side}.FaceAlpha]=fibalpha{:};
                    end
                end

                % Set colorbar tick positions and labels
                if ~isempty(allvals)
                    if domultitract
                        if obj.subscore.special_case
                           if obj.posvisible && obj.negvisible
                                tick{group} = [1, length(fibcmap{group})];
                                poscbvals = sort(allvals(allvals>0));
                                negcbvals = sort(allvals(allvals<0));
                                if ~isempty(negcbvals) && ~isempty(poscbvals)
                                    ticklabel{group} = [negcbvals(1), poscbvals(end)];
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
                                tick{group} = [1, length(fibcmap{group})];
                                poscbvals = sort(allvals(allvals>0));
                                negcbvals = sort(allvals(allvals<0));
                                ticklabel{group} = [negcbvals(1), poscbvals(end)];
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
                            tick{group} = [1, length(fibcmap{group})];
                            poscbvals = sort(allvals(allvals>0));
                            negcbvals = sort(allvals(allvals<0));
                            ticklabel{group} = [negcbvals(1), poscbvals(end)];
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
                if strfind(thisentry,'cleartune')
                    thisentry=strrep(thisentry,'cleartune','');
                    k=strfind(thisentry,'_');
                    ctentry=str2double(thisentry(1:k-1));
                    ctside=str2double(thisentry(k+1:end));
                    weights{ctside}=weights{ctside}+...
                        full(obj.cleartuneresults.(ea_conn2connid(obj.connectome)).(ea_method2methodid(obj)).fibsval{ctside}(:,ctentry));
                elseif strfind(thisentry,'results')
                    weights={ones(size(obj.cleartuneresults.(ea_conn2connid(obj.connectome)).fibcell{1},1),1),...
                        ones(size(obj.cleartuneresults.(ea_conn2connid(obj.connectome)).fibcell{2},1),1)};
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
