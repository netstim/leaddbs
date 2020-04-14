classdef ea_disctract < handle
    % Discriminative fiber class to handle visualizations of discriminative fibers in lead dbs resultfig / 3D Matlab figures
    % A. Horn

    properties (SetObservable)
        M % content of lead group project
        resultfig % figure handle to plot results
        ID % name / ID of discriminative fibers object
        posvisible = 1 % pos tract visible
        negvisible = 0 % neg tract visible
        showposamount = [2 2] % two entries for right and left
        shownegamount = [2 2] % two entries for right and left
        connthreshold = 20
        efieldthreshold = 2500
        statmetric = 1 % entry from discfiber settings as initially specified in prefs.machine.lg (which are separately stored for each analysis/object).
        efieldmetric = 'Peak' % if statmetric == 2, efieldmetric can calculate sum, mean or peak along tracts
        poscolor = [1,0,0] % positive main color
        negcolor = [0,0,1] % negative main color
        splitbygroup = 0
        % the  main content variables:

        results
        % includes subfields by results.connectomename.ttests /
        % results.connectomename.efields with
        % fibcell % cell of all fibers connected, sorted by side
        % fibsval % connection weight value for each fiber to each VTA
        % fibweights % usually T- or R-values associated with each tract
        %
        drawobject % actual streamtube handle
        patientselection % selected patients to include. Note that connected fibers are always sampled from all (& mirrored) VTAs of the lead group file
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
        % stats: (how many fibers available and shown etc for GUI)
        stats
        % additional settings:
        kfold = 5 % divide into k sets when doing k-fold CV
        Nperm = 1000 % how many permutations in leave-nothing-out permtest strategy
        adjustforgroups = 1 % adjust correlations for group effects
    end

    properties (Access = private)
        switchedFromSpace=3 % if switching space, this will protocol where from
    end

    methods
        function obj=ea_disctract() % class constructor

        end

        function initialize(obj,lgpath,resultfig)
            D = load(lgpath);
            if isfield(D, 'M') % Lead Group analysis path loaded
                obj.M = D.M;
                obj.leadgroup = lgpath;

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
            else
                ea_error('You have opened a file of unknown type.')
                return
            end
        end

        function calculate(obj)
            % check that this has not been calculated before:
            if ~isempty(obj.results) % something has been calculated
                if isfield(obj.results,ea_conn2connid(obj.connectome))
                    if isfield(obj.results.(ea_conn2connid(obj.connectome)),ea_method2methodid(obj)) % this combination was already calculated.
                        answ=questdlg('This has already been calculated. Are you sure you want to re-calculate everything?','Recalculate Results','No','Yes','No');
                        if ~strcmp(answ,'Yes')
                            return
                        end
                    end
                end
            end

            cfile=[ea_getconnectomebase('dMRI'),obj.connectome,filesep,'data.mat'];
            [fibcell,fibsin,XYZmm,niivx,valsmm]=ea_discfibers_getfibcell(obj,cfile);
            switch obj.statmetric
                case 1 % ttests
                    [fibsval]=ea_discfibers_heatfibertracts(obj,fibcell,fibsin,XYZmm,niivx);
                    % Main output of results - this is all we will ever need if statmetric
                    % and connectome wont change
                    obj.results.(ea_conn2connid(obj.connectome)).(ea_method2methodid(obj)).fibsval=fibsval;
                case 2 % spearmans R
                    [fibsval_sum,fibsval_mean,fibsval_peak,fibsval_5peak]=ea_discfibers_heatfibertracts_corr(obj,fibcell,XYZmm,niivx,valsmm);
                    obj.results.(ea_conn2connid(obj.connectome)).(ea_method2methodid(obj,'Sum')).fibsval=fibsval_sum;
                    obj.results.(ea_conn2connid(obj.connectome)).(ea_method2methodid(obj,'Mean')).fibsval=fibsval_mean;
                    obj.results.(ea_conn2connid(obj.connectome)).(ea_method2methodid(obj,'Peak')).fibsval=fibsval_peak;
                    obj.results.(ea_conn2connid(obj.connectome)).(ea_method2methodid(obj,'Peak 5%')).fibsval=fibsval_5peak;
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
                    if isempty(obj.M.stats(pt).ea_stats.stimulation.efield(side).volume)
                        val=0;
                    else
                        val=obj.M.stats(pt).ea_stats.stimulation.efield(side).volume;
                    end
                    Efieldmags(pt,side)=val;
                end
            end
        end
        function refreshlg(obj)
            D = load(obj.leadgroup);
            obj.M = D.M;
        end
            
        function coh = getcohortregressor(obj)
            coh=ea_cohortregressor(obj.M.patient.group(obj.patientselection));
        end

        function [I,Ihat] = loocv(obj)
            allpts=obj.patientselection;
            fibsval=obj.results.(ea_conn2connid(obj.connectome)).(ea_method2methodid(obj)).fibsval;
            I=obj.responsevar;
            for side=1:2
                nfibsval{side}=fibsval{side};
                nfibsval{side}(nfibsval{side}==0)=0; % only used in spearmans correlations
            end
            for pt=allpts
                opts=allpts;  % generate variable of patients on which model will be built.
                opts(opts==pt)=[];
                [vals]=ea_disc_calcstats(obj,opts);
                for side=1:2
                    switch obj.statmetric
                        case 1 % ttests, vtas - see Baldermann et al. 2019 Biological Psychiatry
                            Ihat(pt,side)=ea_nansum(vals{1,side}.*fibsval{side}(:,pt)); % I hat is the estimate of improvements (not scaled to real improvements)
                        case 2 % spearmans correlations, efields - see Li et al. 2019 TBP
                            Ihat(pt,side)=ea_nansum(vals{1,side}.*nfibsval{side}(:,pt)); % I hat is the estimate of improvements (not scaled to real improvements)
                    end
                    ea_dispercent(pt/length(allpts));
                end
            end
            if size(obj.responsevar,2)==2 % hemiscores
                Ihat=Ihat(:); % compare hemiscores (electrode wise)
                I=I(:);
            else
                Ihat=ea_nanmean(Ihat,2); % compare bodyscores (patient wise)
            end
            ea_dispercent(1,'end');
            loginx=zeros(size(Ihat));
            loginx(allpts,:)=1;
            Ihat(~loginx)=nan; % make sure info of not included patients are not used
        end

        function [I,Ihat] = lococv(obj)
            [I,Ihat]=groupwisecv(obj,obj.M.patient.group');
        end

        function [I,Ihat] = kfoldcv(obj) 
            % group assignment:
            setlen=floor(length(obj.patientselection)/obj.kfold);
            cnt=1;
            for fold=1:obj.kfold
            groupassignment(cnt:cnt+setlen-1)=fold;
            cnt=cnt+setlen;
            end
            groupassignment(length(groupassignment):(length(groupassignment)+(length(obj.patientselection)-length(groupassignment))))=fold; % fill up remnants
            %randpermute
            groupassignment=groupassignment(randperm(length(groupassignment)));
            [I,Ihat]=groupwisecv(obj,groupassignment);

        end

        function [I,Ihat]=groupwisecv(obj,groupassignment)
            allpts=obj.patientselection;
            fibsval=obj.results.(ea_conn2connid(obj.connectome)).(ea_method2methodid(obj)).fibsval;
            I=obj.responsevar;
            for side=1:2  % only used in spearmans correlations
                if obj.statmetric==2
                    nfibsval{side}=fibsval{side};
                    nfibsval{side}(nfibsval{side}==0)=0;
                end
            end
            if length(unique(groupassignment))<2
               ea_error(['Only 1 set of patients assigned. Cross-validation not possible.']); 
            end
            ea_dispercent(0,'Iterating sets');
            for group=unique(groupassignment)
                opts=allpts; % generate variable of patients on which model will be built.
                opts(groupassignment(allpts)==group)=[];
                [vals]=ea_disc_calcstats(obj,opts);

                for side=1:2
                    switch obj.statmetric
                        case 1 % ttests, vtas - see Baldermann et al. 2019 Biological Psychiatry
                            for pt=find(groupassignment==group)
                                Ihat(pt,side)=ea_nansum(vals{1,side}.*double(fibsval{side}(:,pt))); % I hat is the estimate of improvements (not scaled to real improvements)
                            end
                        case 2 % spearmans correlations, efields - see Li et al. 2019 TBP
                            for pt=find(groupassignment==group)
                                Ihat(pt,side)=ea_nansum(vals{1,side}.*nfibsval{side}(:,pt)); % I hat is the estimate of improvements (not scaled to real improvements)
                            end
                    end
                end
                ea_dispercent(group/length(unique(groupassignment)));
            end
            ea_dispercent(1,'end');
            if size(obj.responsevar,2)==2 % hemiscores
                Ihat=Ihat(:); % compare hemiscores (electrode wise)
                I=I(:);
            else
                Ihat=ea_nanmean(Ihat,2); % compare bodyscores (patient wise)
            end
            ea_dispercent(1,'end');
            loginx=zeros(size(Ihat)); loginx(allpts,:)=1;
            Ihat(~loginx)=nan; % make sure info of not included patients are not used
        end
        
        function [I,Ihat,R0,R1,pperm] = lnopb(obj)
            allpts=obj.patientselection;
            Numperm=obj.Nperm; % run as many as Nperm permutations

            fibsval=obj.results.(ea_conn2connid(obj.connectome)).(ea_method2methodid(obj)).fibsval;
            I=obj.responsevar;
            R0 = zeros(Numperm+1,1);
            for side=1:2
                nfibsval{side}=fibsval{side}; nfibsval{side}(nfibsval{side}==0)=0; % only used in spearmans correlations
            end

            for perm=1:Numperm+1
                Ihat{perm} = zeros(length(I),2);
                if perm>1
                    Iperm=I(randperm(length(I)));
                    fprintf('Calculating permutations: %d/%d\n\n', perm-1, Numperm);
                else % use real empirical set in first run
                    Iperm=I;
                end

                [vals]=ea_disc_calcstats(obj,allpts,Iperm);

                for side=1:2
                    switch obj.statmetric
                        case 1 % ttests, vtas - see Baldermann et al. 2019 Biological Psychiatry
                            for pt=allpts
                                thisptval=fibsval{side}(:,pt); % this patients connections to each fibertract (1 = connected, 0 = unconnected)
                                Ihat{perm}(pt,side)=ea_nansum(vals{side}.*thisptval); % I hat is the estimate of improvements (not scaled to real improvements)
                            end
                        case 2 % spearmans correlations, efields - see Irmen et al. 2019 TBP
                            for pt=allpts
                                Ihat{perm}(pt,side)=ea_nansum(vals{side}.*nfibsval{side}(:,pt)); % I hat is the estimate of improvements (not scaled to real improvements)
                            end
                    end
                    R0(perm,side)=corr(Iperm,Ihat{perm}(:,side),'type','Spearman','rows','pairwise');
                end
            end

            R1=ea_nanmean(R0(1,:),2); % real correlation value when using empirical values
            R0=ea_nanmean(R0(2:end,:),2); % 1-by-Nperm set of R values

            % generate null distribution
            R0=abs(R0);
            R0=sort(R0,'descend');
            p95=R0(round(0.05*Numperm));
            v=ea_searchclosest(R0,R1);

            pperm=v/Numperm;
            disp(['Permuted p = ',sprintf('%0.2f',pperm),'.']);

            Ihat=Ihat{1};
            Ihat=ea_nanmean(Ihat,2);
            loginx=zeros(size(Ihat)); loginx(allpts,:)=1;
            Ihat(~loginx)=nan; % make sure info of not included patients are not used
        end

        function save(obj)
            tractset=obj;
            [pth,fn]=fileparts(tractset.leadgroup);
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
            [vals]=ea_disc_calcstats(obj);
            obj.stats.pos.available(1)=sum(vals{1,1}>0);
            obj.stats.neg.available(1)=sum(vals{1,1}<0);            
            obj.stats.pos.available(2)=sum(vals{1,2}>0);
            obj.stats.neg.available(2)=sum(vals{1,2}<0);

            set(0,'CurrentFigure',obj.resultfig);

            dogroups=size(vals,1)>1; % if color by groups is set will be positive.
            linecols=obj.M.groups.color;
            if isempty(obj.drawobject) % check if prior object has been stored
                obj.drawobject=getappdata(obj.resultfig,['dt_',obj.ID]); % store handle of tract to figure.
            end
            for tract=1:length(obj.drawobject)
                delete(obj.drawobject{tract});

            end
            % reset colorbar
            obj.colorbar=[];
            if ~any([obj.posvisible,obj.negvisible])
                return
            end

            for group=1:size(vals,1) % vals will have 1x2 in case of bipolar drawing and Nx2 in case of group-based drawings (where only positives are shown).
                % Contruct default blue to red colormap
                colormap(gray);
                if dogroups
                    fibcmap{group} = ea_colorgradient(1024, [1,1,1], linecols(group,:));
                    setappdata(obj.resultfig, ['fibcmap',obj.ID], fibcmap);
                else
                    fibcmap{1} = ea_colorgradient(1024, obj.negcolor, [1,1,1], obj.poscolor);
                    setappdata(obj.resultfig, ['fibcmap',obj.ID], fibcmap);
                end

                % Set alphas of fibers with light color to 0
                colorbarThreshold = 0.60; % Percentage of the pos/neg color to be kept
                negUpperBound=ceil(size(fibcmap{group},1)/2*colorbarThreshold);
                poslowerBound=floor((size(fibcmap{group},1)-size(fibcmap{group},1)/2*colorbarThreshold));
                for side=1:2
                    fibcell{group,side}=obj.results.(ea_conn2connid(obj.connectome)).fibcell{side}(~isnan(vals{group,side}));
                    if dogroups % introduce small jitter for visualization
                        fibcell{group,side}=ea_disc_addjitter(fibcell{group,side},0.01);
                    end
                    vals{group,side}=vals{group,side}(~isnan(vals{group,side}))'; % final weights for surviving fibers

                    posvals{group,side}=sort(vals{group,side}(vals{group,side}>0),'descend');
                    negvals{group,side}=sort(vals{group,side}(vals{group,side}<0),'ascend');

                    try
                        posthresh{group,side}=posvals{group,side}(ceil(((obj.showposamount(side)+eps)/100)*length(posvals{group,side})));
                    catch
                        posthresh{group,side}=inf;
                    end
                    try
                        negthresh{group,side}=negvals{group,side}(ceil(((obj.shownegamount(side)+eps)/100)*length(negvals{group,side})));
                    catch
                        negthresh{group,side}=-inf;
                    end
                    % Remove vals and fibers outside the thresholding range
                    remove=logical(logical(vals{group,side}<posthresh{group,side}) .* logical(vals{group,side}>negthresh{group,side}));
                    vals{group,side}(remove)=[];
                    fibcell{group,side}(remove)=[];

                    tvalsRescale{group,side} = vals{group,side};
                    tvalsRescale{group,side}(vals{group,side}>0) = ea_rescale(vals{group,side}(vals{group,side}>0), [0 1]);
                    tvalsRescale{group,side}(vals{group,side}<0) = ea_rescale(vals{group,side}(vals{group,side}<0), [-1 0]);

                    fibcolorInd{group,side}=tvalsRescale{group,side}*(size(fibcmap{group},1)/2-0.5);
                    fibcolorInd{group,side}=fibcolorInd{group,side}+(size(fibcmap{group},1)/2+0.5);

                    alphas{group,side}=zeros(size(fibcolorInd{group,side},2),1);

                    alphas{group,side}(round(fibcolorInd{group,side})>=poslowerBound) = obj.posvisible;
                    alphas{group,side}(round(fibcolorInd{group,side})<=negUpperBound) = obj.negvisible;

                    fibalpha=mat2cell(alphas{group,side},ones(size(fibcolorInd{group,side},2),1));

                    % Plot fibers if any survived
                    if ~isempty(fibcell{group,side})
                        obj.drawobject{group,side}=streamtube(fibcell{group,side},0.2);

                        nones=repmat({'none'},size(fibcolorInd{group,side}));
                        [obj.drawobject{group,side}.EdgeColor]=nones{:};

                        % Calulate fiber colors
                        colors=fibcmap{group}(round(fibcolorInd{group,side}),:);
                        fibcolor=mat2cell(colors,ones(size(fibcolorInd{group,side})));

                        % Set fiber colors and alphas
                        [obj.drawobject{group,side}.FaceColor]=fibcolor{:};
                        [obj.drawobject{group,side}.FaceAlpha]=fibalpha{:};
                    end
                end
                
                obj.stats.pos.shown(1)=sum(vals{1,1}>0);
                obj.stats.neg.shown(1)=sum(vals{1,1}<0);
                obj.stats.pos.shown(2)=sum(vals{1,2}>0);
                obj.stats.neg.shown(2)=sum(vals{1,2}<0);

                % Set colorbar tick positions and labels
                if ~any([isempty(vals{group,1}),isempty(vals{group,2})])
                    cbvals = [vals{group,1}(logical(alphas{group,1})),vals{group,2}(logical(alphas{group,2}))];
                    % cbvals=tvalsRescale{group,side}(logical(alphas));
                    if obj.posvisible && ~obj.negvisible
                        cbmap{group} = fibcmap{group}(ceil(length(fibcmap{group})/2+0.5):end,:);
                        tick{group} = [poslowerBound, length(fibcmap{group})] - floor(length(fibcmap{group})/2) ;
                        poscbvals = sort(cbvals(cbvals>0));
                        ticklabel{group} = [poscbvals(1), poscbvals(end)];
                        ticklabel{group} = arrayfun(@(x) num2str(x,'%.2f'), ticklabel{group}, 'Uni', 0);
                    elseif ~obj.posvisible && obj.negvisible
                        cbmap{group} = fibcmap{group}(1:floor(length(fibcmap{group})/2-0.5),:);
                        tick{group} = [1, negUpperBound];
                        negcbvals = sort(cbvals(cbvals<0));
                        ticklabel{group} = [negcbvals(1), negcbvals(end)];
                        ticklabel{group} = arrayfun(@(x) num2str(x,'%.2f'), ticklabel{group}, 'Uni', 0);
                    elseif obj.posvisible && obj.negvisible
                        cbmap{group} = fibcmap{group};
                        tick{group} = [1, negUpperBound, poslowerBound, length(fibcmap{group})];
                        poscbvals = sort(cbvals(cbvals>0));
                        negcbvals = sort(cbvals(cbvals<0));
                        ticklabel{group} = [min(cbvals), negcbvals(end), poscbvals(1), max(cbvals)];
                        ticklabel{group} = arrayfun(@(x) num2str(x,'%.2f'), ticklabel{group}, 'Uni', 0);
                    end
                end
            end

            % store colorbar in object
            if exist('cbmap','var') % could be no fibers present at all.
                obj.colorbar.cmap = cbmap;
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
