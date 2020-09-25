classdef ea_disctract < handle
    % Discriminative fiber class to handle visualizations of discriminative fibers in lead dbs resultfig / 3D Matlab figures
    % A. Horn

    properties (SetObservable)
        M % content of lead group project
        resultfig % figure handle to plot results
        ID % name / ID of discriminative fibers object
        posvisible = 1 % pos tract visible
        negvisible = 0 % neg tract visible
        showposamount = [25 25] % two entries for right and left
        shownegamount = [25 25] % two entries for right and left
        connthreshold = 20
        efieldthreshold = 2500
        statmetric = 1 % stats metric to use, 1 = ttest, 2 = spearman
        efieldmetric = 'Peak' % if statmetric == 2, efieldmetric can calculate sum, mean or peak along tracts
        poscolor = [1,0,0] % positive main color
        negcolor = [0,0,1] % negative main color
        splitbygroup = 0

        results
        % Subfields:
        % results.(connectomename).fibcell: cell of all fibers connected, sorted by side
        % results.(connectomename).ttests.fibsval % connection status for each fiber to each VTA
        % results.(connectomename).spearman_sum.fibsval % connection weights for each fiber to each VTA
        % results.(connectomename).spearman_mean.fibsval % connection weights for each fiber to each VTA
        % results.(connectomename).spearman_peak.fibsval % connection weights for each fiber to each VTA
        % results.(connectomename).spearman_5peak.fibsval % connection weights for each fiber to each VTA

        fiberdrawn % struct contains fibercell and vals drawn in the resultfig
        drawobject % actual streamtube handle
        patientselection % selected patients to include. Note that connected fibers are always sampled from all (& mirrored) VTAs of the lead group file
        customselection % selected patients in the custom test list
        allpatients % list of all patients (as from M.patient.list)
        mirrorsides = 0 % flag to mirror VTAs / Efields to contralateral sides using ea_flip_lr_nonlinear()
        responsevar % response variable
        responsevarlabel % label of response variable
        multresponsevarneg = 0 % multiply response variable by -1
        covars = {} % covariates
        covarlabels = {} % covariate labels
        analysispath % where to store results
        leadgroup % redundancy protocol only, path to original lead group project
        connectome % redundancy protocol only, name of underlying connectome
        colorbar % colorbar information
        % stats: (how many fibers available and shown etc for GUI)
        stats
        % additional settings:
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
            D = load(datapath);
            if isfield(D, 'M') % Lead Group analysis path loaded
                obj.M = D.M;
                obj.leadgroup = datapath;

                testID = obj.M.guid;
                ea_mkdir([fileparts(obj.leadgroup),filesep,'disctracts',filesep]);
                id = 1;
                while exist([fileparts(obj.leadgroup),filesep,'disctracts',filesep,testID,'.mat'],'file')
                    testID = [obj.M.guid, '_', num2str(id)];
                    id = id + 1;
                end
                obj.ID = testID;
                obj.resultfig = resultfig;
                obj.allpatients = obj.M.patient.list;
                obj.patientselection = obj.M.ui.listselect;
                obj.responsevarlabel = obj.M.clinical.labels{1};
                obj.covarlabels={'Stimulation Amplitude'};
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

            cfile = [ea_getconnectomebase('dMRI'), obj.connectome, filesep, 'data.mat'];
            vatlist = ea_discfibers_getvats(obj);
            [fibsvalBin, fibsvalSum, fibsvalMean, fibsvalPeak, fibsval5Peak, fibcell] = ea_discfibers_calcvals(vatlist, cfile);

            obj.results.(ea_conn2connid(obj.connectome)).('ttests').fibsval = fibsvalBin;
            obj.results.(ea_conn2connid(obj.connectome)).('spearman_sum').fibsval = fibsvalSum;
            obj.results.(ea_conn2connid(obj.connectome)).('spearman_mean').fibsval = fibsvalMean;
            obj.results.(ea_conn2connid(obj.connectome)).('spearman_peak').fibsval = fibsvalPeak;
            obj.results.(ea_conn2connid(obj.connectome)).('spearman_5peak').fibsval = fibsval5Peak;
            obj.results.(ea_conn2connid(obj.connectome)).fibcell = fibcell;
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
               msgbox('LEAD_groupanalysis file has vanished. Please select file.');
               [fn,pth]=uigetfile();
               obj.leadgroup=fullfile(pth,fn);
            end
            D = load(obj.leadgroup);
            obj.M = D.M;
        end

        function coh = getcohortregressor(obj)
            coh=ea_cohortregressor(obj.M.patient.group(obj.patientselection));
        end

        function [I, Ihat] = loocv(obj)
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
            cvp = cvpartition(length(obj.patientselection), 'KFold', obj.kfold);
            [I, Ihat] = crossval(obj, cvp);
        end

        function [I, Ihat] = lno(obj, Iperm)
            cvp = cvpartition(length(obj.patientselection), 'resubstitution');
            if ~exist('Iperm', 'var')
                [I, Ihat] = crossval(obj, cvp);
            else
                [I, Ihat] = crossval(obj, cvp, Iperm);
            end
        end

        function [I, Ihat] = crossval(obj, cvp, Iperm)
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
                I = obj.responsevar(patientsel);
            else
                I = Iperm(patientsel);
            end

            % Ihat is the estimate of improvements (not scaled to real improvements)
            Ihat = nan(length(patientsel),2);

            fibsval = full(obj.results.(ea_conn2connid(obj.connectome)).(ea_method2methodid(obj)).fibsval);
            for side=1:2  % only used in spearmans correlations
                if obj.statmetric==2
                    nfibsval{side} = fibsval{side};
                    nfibsval{side}(nfibsval{side}==0) = 0;
                end
            end

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
                    [vals,~,usedidx] = ea_discfibers_calcstats(obj, patientsel(training));
                else
                    [vals,~,usedidx] = ea_discfibers_calcstats(obj, patientsel(training), Iperm);
                end

                for side=1:2
                    if ~isempty(vals{1,side})
                        switch obj.statmetric
                            case 1 % ttests, vtas - see Baldermann et al. 2019 Biological Psychiatry
                            	Ihat(test,side) = ea_nansum(vals{1,side}.*fibsval{1,side}(usedidx{1,side},patientsel(test)));
                            case 2 % spearmans correlations, efields - see Irmen et al. 2020 Annals of Neurology
                            	Ihat(test,side) = ea_nansum(vals{1,side}.*nfibsval{side}(usedidx{1,side},patientsel(test)));
                        end
                    end
                end
            end

            if cvp.NumTestSets == 1
                Ihat = Ihat(test,:);
                I = I(test);
            end

            if size(obj.responsevar,2)==2 % hemiscores
                Ihat = Ihat(:); % compare hemiscores (electrode wise)
                I = I(:);
            else
                Ihat = ea_nanmean(Ihat,2); % compare bodyscores (patient wise)
            end
        end

        function [Iperm, Ihat, R0, R1, pperm, Rp95] = lnopb(obj, corrType)
            if ~exist('corrType', 'var')
                corrType = 'Spearman';
            end

            fibsval = full(obj.results.(ea_conn2connid(obj.connectome)).(ea_method2methodid(obj)).fibsval);
            for side=1:2
                nfibsval{side}=fibsval{side};
                nfibsval{side}(nfibsval{side}==0)=0; % only used in spearmans correlations
            end

            numPerm = obj.Nperm;

            Iperm = ea_shuffle(obj.responsevar, numPerm, obj.patientselection)';
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
            R0 = sort(abs(R(2:end)),'descend');
            Rp95 = R0(round(0.05*numPerm));
            v = ea_searchclosest(R0, R1);
            pperm = v/numPerm;
            disp(['Permuted p = ',sprintf('%0.2f',pperm),'.']);

            % Return only selected I
            Iperm = Iperm(obj.patientselection,:);
        end

        function save(obj)
            tractset=obj;
            pth = fileparts(tractset.leadgroup);
            tractset.analysispath=[pth,filesep,'disctracts',filesep,obj.ID,'.mat'];
            ea_mkdir([pth,filesep,'disctracts']);
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

        function draw(obj)
            [vals,fibcell]=ea_discfibers_calcstats(obj);
            obj.fiberdrawn.fibcell = fibcell;
            obj.fiberdrawn.vals = vals;

            obj.stats.pos.shown(1)=sum(vals{1,1}>0);
            obj.stats.neg.shown(1)=sum(vals{1,1}<0);
            obj.stats.pos.shown(2)=sum(vals{1,2}>0);
            obj.stats.neg.shown(2)=sum(vals{1,2}<0);

            set(0,'CurrentFigure',obj.resultfig);

            dogroups=size(vals,1)>1; % if color by groups is set will be positive.
            linecols=obj.M.groups.color;
            if isempty(obj.drawobject) % check if prior object has been stored
                obj.drawobject=getappdata(obj.resultfig,['dt_',obj.ID]); % store handle of tract to figure.
            end
            for tract=1:numel(obj.drawobject)
                delete(obj.drawobject{tract});
            end

            % reset colorbar
            obj.colorbar=[];
            if ~any([obj.posvisible,obj.negvisible])
                return
            end

            for group=1:size(vals,1) % vals will have 1x2 in case of bipolar drawing and Nx2 in case of group-based drawings (where only positives are shown).
                % Contruct default blue to red colormap
                allvals = vertcat(vals{group,:});
                if isempty(allvals)
                    continue;
                end
                colormap(gray);
                gradientLevel = 1024;
                cmapShiftRatio = 0.5;
                shiftedCmapStart = round(gradientLevel*cmapShiftRatio)+1;
                shiftedCmapEnd = gradientLevel-round(gradientLevel*cmapShiftRatio);
                shiftedCmapLeftEnd = gradientLevel/2-round(gradientLevel/2*cmapShiftRatio);
                shiftedCmapRightStart = round(gradientLevel/2*cmapShiftRatio)+1;
                if dogroups
                    if obj.posvisible && ~obj.negvisible
                        cmap = ea_colorgradient(gradientLevel, [1,1,1], linecols(group,:));
                        fibcmap{group} = ea_colorgradient(gradientLevel, cmap(shiftedCmapStart,:), linecols(group,:));
                        cmapind = round(normalize(allvals,'range',[1,gradientLevel]));
                        alphaind = ones(size(allvals));
                        % alphaind = normalize(allvals, 'range');
                    elseif ~obj.posvisible && obj.negvisible
                        cmap = ea_colorgradient(gradientLevel, linecols(group,:), [1,1,1]);
                        fibcmap{group} = ea_colorgradient(gradientLevel, linecols(group,:), cmap(shiftedCmapEnd,:));
                        cmapind = round(normalize(allvals,'range',[1,gradientLevel]));
                        alphaind = ones(size(allvals));
                        % alphaind = normalize(-allvals, 'range');
                    else
                        warndlg(sprintf(['Please choose either "Show Positive Fibers" or "Show Negative Fibers".',...
                            '\nShow both positive and negative fibers is not supported when "Color by Group Variable" is on.']));
                        return;
                    end
                else
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

                cmapind = mat2cell(cmapind, [numel(vals{group,1}), numel(vals{group,2})])';
                alphaind = mat2cell(alphaind, [numel(vals{group,1}), numel(vals{group,2})])';

                for side=1:2
                    if dogroups % introduce small jitter for visualization
                        fibcell{group,side}=ea_discfibers_addjitter(fibcell{group,side},0.01);
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
    end

    methods (Static)
        function changeevent(~,event)
            update_trajectory(event.AffectedObject,event.Source.Name);
        end
    end
end
