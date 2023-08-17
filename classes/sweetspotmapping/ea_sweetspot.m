classdef ea_sweetspot < handle
    % Sweetspot analysis class to handle visualizations of sweetspots in lead dbs resultfig / 3D Matlab figures
    % A. Horn

    properties (SetObservable)
        M % content of lead group project
        resultfig % figure handle to plot results
        ID % name / ID of sweetspot object
        posvisible = 1 % sweetspot visible
        negvisible = 0 % sourspot visible

        efieldthreshold = 200
        statlevel = 'VTAs' % stats metric to use, 1 = active coordinates, 2 = efields, 3 = vtas
        stattest = 'T-Test';
        stat0hypothesis = 'Zero';
        statimpthreshold = 0;
        statNthreshold = 0;
        statamplitudecorrection = 'None';
        statnormalization = 'None';
        corrtype = 'Spearman' % correlation strategy in case of using E-Fields.
        coverthreshold = 20; % of vtas needed to cover a single voxel to be considered
        posBaseColor = [1,1,1] % positive main color
        posPeakColor = [0.9176,0.2000,0.1373] % positive peak color

        negBaseColor = [1,1,1] % negative main color
        negPeakColor = [0.2824,0.6157,0.9725] % negative peak color

        splitbygroup = 0
        showsignificantonly = 1
        alphalevel = 0.05
        multcompstrategy = 'Uncorrected'; % could be 'Bonferroni'

        autorefresh=1;

        results
        % Subfields:
        cvlivevisualize = 0; % if set to 1 shows crossvalidation results during processing.
        basepredictionon = 'Mean of Scores';
        spotdrawn % struct contains sweetspot drawn in the resultfig
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
        colorbar % colorbar information
        % stats: (how many fibers available and shown etc for GUI)
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
        function obj=ea_sweetspot(analysispath) % class constructor
            if exist('analysispath', 'var') && ~isempty(analysispath)
                obj.analysispath = analysispath;
                [~, ID] = fileparts(obj.analysispath);
                obj.ID = ID;
            end
        end

        function initialize(obj,datapath,resultfig)
            D = load(datapath, '-mat');
            if isfield(D, 'M') % Lead Group analysis path loaded
                obj.M = D.M;
                obj.leadgroup = datapath;

                testID = obj.M.guid;
                ea_mkdir([fileparts(obj.leadgroup),filesep,'sweetspots',filesep]);
                id = 1;
                while exist([fileparts(obj.leadgroup),filesep,'sweetspots',filesep,testID,'.sweetspot'],'file')
                    testID = [obj.M.guid, '_', num2str(id)];
                    id = id + 1;
                end
                obj.ID = testID;
                obj.resultfig = resultfig;

                if isfield(obj.M,'pseudoM')
                    obj.allpatients = obj.M.ROI.list;
                    obj.patientselection = 1:length(obj.M.ROI.list);
                    obj.M.root = [fileparts(datapath),filesep];
                    %obj.M.patient.list=obj.M.ROI.list; % copies
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
                obj.covarlabels={};
            elseif  isfield(D, 'sweetspot')  % Saved sweetspot class loaded
                props = properties(D.sweetspot);
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
            obj.calculate;
        end

        function calculate(obj)
            % in case of the sweetspot explorer, calculate rather means to
            % gather all E-Fields. To keep consistency of the logic with
            % discfiberexplorer and networkmappingexplorer, we will keep
            % the same name (calculate) for the function, nonetheless.

            % check that results aren't already there
            if ~isempty(obj.results) % vtas already gathered in
                return
            end

            if isfield(obj.M,'pseudoM')
                vatlist = obj.M.ROI.list;
            else
                vatlist = ea_sweetspot_getvats(obj);
            end
            [AllX,space] = ea_exportefieldmapping(vatlist,obj);

            obj.results.efield = AllX;
            obj.results.space = space;

            if ~isfield(obj.M,'pseudoM')
                % get active coordinates, as well
                for pt=1:length(obj.M.patient.list)
                    for side=1:2
                        obj.results.activecnt{side}(pt,:)=...
                            mean(obj.M.elstruct(pt).coords_mm{side}(find(obj.M.S(pt).activecontacts{side}),:),1); %#ok<FNDSB> % find is necessary here
                    end
                end
                for side=1:2
                    obj.results.activecnt{side}=[obj.results.activecnt{side};ea_flip_lr_nonlinear(obj.results.activecnt{side})];
                end
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
                I = obj.responsevar(patientsel,:);
            else
                I = Iperm(patientsel,:);
            end

            % Ihat is the estimate of improvements (not scaled to real improvements)
            Ihat = nan(length(patientsel),2);

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
                        [vals] = ea_sweetspot_calcstats(obj, patientsel(training));
                        obj.draw(vals);
                        drawnow;
                    else
                        [vals] = ea_sweetspot_calcstats(obj, patientsel(training));
                    end
                else
                    if obj.cvlivevisualize
                        [vals] = ea_sweetspot_calcstats(obj, patientsel(training), Iperm);
                        obj.draw(vals);
                        drawnow;
                    else
                        [vals] = ea_sweetspot_calcstats(obj, patientsel(training), Iperm);
                    end
                end

                for side=1:numel(vals)
                    if ~isempty(vals{1,side})
                        switch obj.statlevel % also differentiate between methods in the prediction part.
                            case 'VTAs'
                                efield = obj.results.efield{side}(patientsel(test),:)';
                                efield(~isnan(efield)) = efield(~isnan(efield)) > obj.efieldthreshold;
                                switch lower(obj.basepredictionon)
                                    case 'mean of scores'
                                        Ihat(test,side) = ea_nanmean(vals{1,side}.*efield,1);
                                    case 'sum of scores'
                                        Ihat(test,side) = ea_nansum(vals{1,side}.*efield,1);
                                    case 'peak of scores'
                                        Ihat(test,side) = ea_discfibers_getpeak(vals{1,side}.*efield, obj.posvisible, obj.negvisible, 'peak');
                                    case 'peak 5% of scores'
                                        Ihat(test,side) = ea_discfibers_getpeak(vals{1,side}.*efield, obj.posvisible, obj.negvisible, 'peak5');
                                end
                            case 'E-Fields'
                                switch lower(obj.basepredictionon)
                                    case 'profile of scores: spearman'
                                        Ihat(test,side) = atanh(ea_corr(vals{1,side},obj.results.efield{side}(patientsel(test),:)','spearman'));
                                    case 'profile of scores: pearson'
                                        Ihat(test,side) = atanh(ea_corr(vals{1,side},obj.results.efield{side}(patientsel(test),:)','pearson'));
                                   case 'profile of scores: bend'
                                        Ihat(test,side) = atanh(ea_corr(vals{1,side},obj.results.efield{side}(patientsel(test),:)','bend'));
                                    case 'mean of scores'
                                        Ihat(test,side) = ea_nanmean(vals{1,side}.*obj.results.efield{side}(patientsel(test),:)',1);
                                    case 'sum of scores'
                                        Ihat(test,side) = ea_nansum(vals{1,side}.*obj.results.efield{side}(patientsel(test),:)',1);
                                    case 'peak of scores'
                                        Ihat(test,side) = ea_discfibers_getpeak(vals{1,side}.*obj.results.efield{side}(patientsel(test),:)', obj.posvisible, obj.negvisible, 'peak');
                                    case 'peak 5% of scores'
                                        Ihat(test,side) = ea_discfibers_getpeak(vals{1,side}.*obj.results.efield{side}(patientsel(test),:)', obj.posvisible, obj.negvisible, 'peak5');
                                end
                        end
                    end
                end
            end

            % restore original view in case of live drawing
            if obj.cvlivevisualize
                obj.draw;
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
            R0 = sort(abs(R(2:end)),'descend');
            Rp95 = R0(round(0.05*numPerm));
            pperm = mean(abs(R0)>=abs(R1));
            disp(['Permuted p = ',sprintf('%0.2f',pperm),'.']);

            % Return only selected I
            Iperm = Iperm(obj.patientselection,:);
        end

        function save(obj)
            sweetspot=obj;
            pth = fileparts(sweetspot.leadgroup);
            sweetspot.analysispath=[pth,filesep,'sweetspots',filesep,obj.ID,'.sweetspot'];
            ea_mkdir([pth,filesep,'sweetspots']);
            rf=obj.resultfig; % need to stash fig handle for saving.
            rd=obj.drawobject; % need to stash handle of drawing before saving.
            try % could be figure is already closed.
                setappdata(rf,['dt_',sweetspot.ID],rd); % store handle of tract to figure.
            end
            sweetspot.resultfig=[]; % rm figure handle before saving.
            sweetspot.drawobject=[]; % rm drawobject.
            save(sweetspot.analysispath,'sweetspot','-v7.3');
            obj.resultfig=rf;
            obj.drawobject=rd;
        end

        function export=draw(obj,vals)
            if ~exist('vals','var')
                [vals]=ea_sweetspot_calcstats(obj);
            end
            obj.spotdrawn.vals=vals;

            obj.stats.pos.shown(1)=sum(vals{1,1}>0);
            obj.stats.neg.shown(1)=sum(vals{1,1}<0);

            set(0,'CurrentFigure',obj.resultfig);

            dogroups=size(vals,1)>1; % if color by groups is set will be positive.
            if ~isfield(obj.M,'groups')
                obj.M.groups.group=ones(length(obj.M.patient.list),1);
                obj.M.groups.color=ea_color_wes('all');
            end
            linecols=obj.M.groups.color;
            if isempty(obj.drawobject) % check if prior object has been stored
                obj.drawobject=getappdata(obj.resultfig,['dt_',obj.ID]); % store handle of tract to figure.
            end
            for s=1:numel(obj.drawobject)
                for ins=1:numel(obj.drawobject{s})
                    try delete(obj.drawobject{s}{ins}.toggleH); end
                    try delete(obj.drawobject{s}{ins}.patchH); end
                    try delete(obj.drawobject{s}{ins}); end
                end
            end
            obj.drawobject={};

            % reset colorbar
            obj.colorbar=[];
            if ~any([obj.posvisible,obj.negvisible])
                export=nan;
            end

            for group=1:size(vals,1) % vals will have 1x2 in case of bipolar drawing and Nx2 in case of group-based drawings (where only positives are shown).
                % Vertcat all values for colorbar construction
                allvals = vertcat(vals{group,:});
                if isempty(allvals) || all(isnan(allvals))
                    ea_cprintf('CmdWinWarnings', 'Empty or all-nan value found!\n');
                    continue;
                else
                    allvals(isnan(allvals)) = 0;
                end

                if obj.posvisible && all(allvals<=0)
                    obj.posvisible = 0;
                    fprintf('\n')
                    warning('off', 'backtrace');
                    warning('No positive values found, posvisible is set to 0 now!');
                    warning('on', 'backtrace');
                    fprintf('\n')
                end

                if obj.negvisible && all(allvals>=0)
                    obj.negvisible = 0;
                    fprintf('\n')
                    warning('off', 'backtrace');
                    warning('No negative values found, negvisible is set to 0 now!');
                    warning('on', 'backtrace');
                    fprintf('\n')
                end

                colormap(gray);
                gradientLevel = length(gray);

                if dogroups
                    if obj.posvisible && ~obj.negvisible
                        voxcmap{group} = ea_colorgradient(gradientLevel, obj.posBaseColor, linecols(group,:));
                    elseif ~obj.posvisible && obj.negvisible
                        voxcmap{group} = ea_colorgradient(gradientLevel, linecols(group,:), obj.negBaseColor);
                    else
                        warndlg(sprintf(['Please choose either "Show Positive Regions" or "Show Negative Regions".',...
                            '\nShow both positive and negative regions is not supported when "Color by Group Variable" is on.']));
                        return;
                    end
                else
                    if obj.posvisible && obj.negvisible
                        cmapLeft = ea_colorgradient(gradientLevel/2, obj.negPeakColor, obj.negBaseColor);
                        cmapRight = ea_colorgradient(gradientLevel/2, obj.posBaseColor, obj.posPeakColor);
                        voxcmap{group} = [cmapLeft;cmapRight];
                    elseif obj.posvisible
                        voxcmap{group} = ea_colorgradient(gradientLevel, obj.posBaseColor, obj.posPeakColor);
                    elseif obj.negvisible
                        voxcmap{group} = ea_colorgradient(gradientLevel, obj.negPeakColor, obj.negBaseColor);
                    end
                end

                for side=1:size(vals,2)
                    res=obj.results.space{side};
                    res.dt(1) = 16;
                    res.img(:)=nan;
                    % Plot voxels if any survived
                    if obj.posvisible
                        % plot positives:
                        posvox=res;
                        posvox.img(:)=0;
                        posvox.img(vals{group,side}>0)=vals{group,side}(vals{group,side}>0);

                        pobj.nii=posvox;
                        pobj.name='Positive';
                        pobj.niftiFilename='Positive.nii';
                        pobj.binary=0;
                        pobj.usesolidcolor=0;
                        pobj.color=obj.posPeakColor;
                        pobj.colormap=ea_colorgradient(gradientLevel, obj.posBaseColor, obj.posPeakColor);
                        pobj.smooth=10;
                        pobj.hullsimplify=0.5;
                        pobj.threshold=0;
                        obj.drawobject{group,side}{1}=ea_roi('Positive.nii',pobj);

                        res=posvox; % keep copy for export
                    end

                    if obj.negvisible
                        % plot negatives:
                        negvox=res;
                        negvox.img(:)=0;
                        negvox.img(vals{group,side}<0)=-vals{group,side}(vals{group,side}<0);

                        pobj.nii=negvox;
                        pobj.name='Negative';
                        pobj.niftiFilename='Negative.nii';
                        pobj.binary=0;
                        pobj.usesolidcolor=0;
                        pobj.color=obj.negPeakColor;
                        pobj.colormap=ea_colorgradient(gradientLevel, obj.negPeakColor, obj.negBaseColor);
                        pobj.smooth=10;
                        pobj.hullsimplify=0.5;
                        pobj.threshold=0;
                        obj.drawobject{group,side}{2}=ea_roi('Negative.nii',pobj);

                        res.img(:)=nansum([res.img(:),-negvox.img(:)],2); % keep copy for export.
                    end
                    res.img(res.img==0)=nan;
                    export{side}=res;
                end

                % Set colorbar tick positions and labels
                if ~isempty(allvals)
                    if obj.posvisible && obj.negvisible
                        tick{group} = [1, gradientLevel/2-10, gradientLevel/2+11, length(voxcmap{group})];
                        poscbvals = sort(allvals(allvals>0));
                        negcbvals = sort(allvals(allvals<0));
                        ticklabel{group} = [negcbvals(1), negcbvals(end), poscbvals(1), poscbvals(end)];
                        ticklabel{group} = arrayfun(@(x) num2str(x,'%.2f'), ticklabel{group}, 'Uni', 0);
                    elseif obj.posvisible
                        tick{group} = [1, length(voxcmap{group})];
                        posvals = sort(allvals(allvals>0));
                        ticklabel{group} = [posvals(1), posvals(end)];
                        ticklabel{group} = arrayfun(@(x) num2str(x,'%.2f'), ticklabel{group}, 'Uni', 0);
                    elseif obj.negvisible
                        tick{group} = [1, length(voxcmap{group})];
                        negvals = sort(allvals(allvals<0));
                        ticklabel{group} = [negvals(1), negvals(end)];
                        ticklabel{group} = arrayfun(@(x) num2str(x,'%.2f'), ticklabel{group}, 'Uni', 0);
                    end
                end
            end

            if ~exist('export','var') % all empty
                for side=1:size(vals,2)
                    res=obj.results.space{side};
                    res.dt(1) = 16;
                    res.img(:)=nan;
                    export{side}=res;
                end
            end

            setappdata(obj.resultfig,['dt_',obj.ID],obj.drawobject); % store handle of surf to figure.

            % store colorbar in object
            if exist('voxcmap','var')
                setappdata(obj.resultfig, ['voxcmap',obj.ID], voxcmap);
                obj.colorbar.cmap = voxcmap;
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
