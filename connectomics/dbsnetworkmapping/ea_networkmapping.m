classdef ea_networkmapping < handle
    % DBS network mapping class to handle visualizations of DBS network mapping analyses in lead dbs resultfig / 3D Matlab figures
    % A. Horn
    
    properties (SetObservable)
        M % content of lead group project
        resultfig % figure handle to plot results
        ID % name / ID of discriminative fibers object
        posvisible = 1 % pos voxels visible
        negvisible = 0 % neg voxels visible
        showposamount = [25 25] % two entries for right and left
        shownegamount = [25 25] % two entries for right and left
        statmetric = 'Correlations (R-map)' % Statistical model to use
        corrtype = 'Spearman' % correlation strategy in case of statmetric == 2.
        poscolor = [0.99,0.75,0.06] % positive main color
        negcolor = [0.15,0.77,0.95] % negative main color
        showsignificantonly = 0
        alphalevel = 0.05
        multcompstrategy = 'FDR'; % could be 'Bonferroni'
        splitbygroup = 0;
        vizmode='Regions'; % way to plot results
        model='Smoothed'; % in case of surface above, on which surface to plot.
        modelLH=1; % show left hemisphere
        modelRH=1; % show right hemisphere
        
        results % stores fingerprints
        outputspace = '222';
        cvlivevisualize = 0; % if set to 1 shows crossvalidation results during processing.
        basepredictionon = 'Spearman Correlations';
        surfdrawn % struct contains data that is actually drawn.
        drawobject % handle to surface drawn on figure
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
        cvmask = 'Gray Matter';
        analysispath % where to store results
        leadgroup % redundancy protocol only, path to original lead group project
        connectome % redundancy protocol only, name of underlying connectome
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
            elseif  isfield(D, 'networkmapping')  % Saved networkmapping class loaded
                props = properties(D.networkmapping);
                for p =  1:length(props) %copy all public properties
                    if ~(strcmp(props{p}, 'analysispath') && ~isempty(obj.analysispath) ...
                            || strcmp(props{p}, 'ID') && ~isempty(obj.ID))
                        obj.(props{p}) = D.networkmapping.(props{p});
                        
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
            
            
            vatlist = ea_networkmapping_getvats(obj);
            [AllX] = ea_networkmapping_calcvals(vatlist, obj.connectome);
            
            obj.results.(ea_conn2connid(obj.connectome)).connval = AllX;
            
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
                I = obj.responsevar(patientsel);
            else
                I = Iperm(patientsel);
            end
            
            % Ihat is the estimate of improvements (not scaled to real improvements)
            Ihat = nan(length(patientsel),1);
            
            connval = full(obj.results.(ea_conn2connid(obj.connectome)).connval);
            
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
                        [vals] = ea_networkmapping_calcstats(obj, patientsel(training));
                        obj.draw(vals);
                        drawnow;
                    else
                        [vals] = ea_networkmapping_calcstats(obj, patientsel(training));
                    end
                else
                    if obj.cvlivevisualize
                        [vals] = ea_networkmapping_calcstats(obj, patientsel(training), Iperm);
                        obj.draw(vals);
                        drawnow;
                    else
                        [vals] = ea_networkmapping_calcstats(obj, patientsel(training), Iperm);
                    end
                end
                
                switch lower(obj.basepredictionon)
                    case 'spearman correlations'
                        Ihat(test) = corr(vals{1}(ea_getmask(ea_mask2maskn(obj)))',...
                            connval(patientsel(test),ea_getmask(ea_mask2maskn(obj)))','rows','pairwise','type','Spearman');
                    case 'pearson correlations'
                        Ihat(test) = corr(vals{1}(ea_getmask(ea_mask2maskn(obj)))',...
                            connval(patientsel(test),ea_getmask(ea_mask2maskn(obj)))','rows','pairwise','type','Pearson');
                    case 'bend correlations'
                        Ihat(test) = ea_bendcorr(vals{1}(ea_getmask(ea_mask2maskn(obj)))',...
                            connval(patientsel(test),ea_getmask(ea_mask2maskn(obj)))');
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
        
        function maskn=ea_mask2maskn(obj)
            switch obj.cvmask
                case 'Gray Matter'
                    switch obj.outputspace
                        case '222'
                            maskn='gray';
                        case '111'
                            maskn='gray_hd';
                        case '555'
                            maskn='gray_5';
                    end
                case 'Brain'
                    switch obj.outputspace
                        case '222'
                            maskn='brain';
                        case '111'
                            maskn='brain_hd';
                        case '555'
                            maskn='brain_5';
                    end
                case 'Cortex & Cerebellum'
                    switch obj.outputspace
                        case '222'
                            maskn='cortexcb';
                        case '111'
                            maskn='cortexcb_hd';
                        case '555'
                            ea_error('Cortex & Cerebellum Mask not supported for 0.5 mm resolution space');
                    end
                case 'Cortex'
                    switch obj.outputspace
                        case '222'
                            maskn='cortex';
                        case '111'
                            maskn='cortex_hd';
                        case '555'
                            maskn='cortex_5';
                    end
                case 'Cerebellum'
                    switch obj.outputspace
                        case '222'
                            maskn='cb';
                        case '111'
                            maskn='cb_hd';
                        case '555'
                            ea_error('Cerebellum Mask not supported for 0.5 mm resolution space');
                    end
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
            v = ea_searchclosest(R0, R1);
            pperm = v/numPerm;
            disp(['Permuted p = ',sprintf('%0.2f',pperm),'.']);
            
            % Return only selected I
            Iperm = Iperm(obj.patientselection,:);
        end
        
        function save(obj)
            networkmapping=obj;
            pth = fileparts(networkmapping.leadgroup);
            networkmapping.analysispath=[pth,filesep,'networkmapping',filesep,obj.ID,'.mat'];
            ea_mkdir([pth,filesep,'networkmapping']);
            rf=obj.resultfig; % need to stash fig handle for saving.
            try % could be figure is already closed.
                setappdata(rf,['dt_',networkmapping.ID],rd); % store handle of tract to figure.
            end
            networkmapping.resultfig=[]; % rm figure handle before saving.
            save(networkmapping.analysispath,'networkmapping','-v7.3');
            obj.resultfig=rf;
        end
        
        
        function draw(obj,vals)
            if ~exist('vals','var')
                [vals]=ea_networkmapping_calcstats(obj);
            end
            obj.surfdrawn.vals=vals;
            
            obj.stats.pos.shown(1)=sum(vals{1,1}>0);
            obj.stats.neg.shown(1)=sum(vals{1,1}<0);
            
            
            set(0,'CurrentFigure',obj.resultfig);
            
            dogroups=size(vals,1)>1; % if color by groups is set will be positive.
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
                return
            end
            
            space=ea_load_nii([ea_getearoot,'templates',filesep,'spacedefinitions',filesep,obj.outputspace,'.nii.gz']);
            for group=1:size(vals,1) % vals will have 1x2 in case of bipolar drawing and Nx2 in case of group-based drawings (where only positives are shown).
                % Contruct default blue to red colormap
                allvals = horzcat(vals{group,:})';
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
                        voxcmap{group} = ea_colorgradient(gradientLevel, cmap(shiftedCmapStart,:), linecols(group,:));
                        cmapind = round(normalize(allvals,'range',[1,gradientLevel]));
                        alphaind = ones(size(allvals));
                        % alphaind = normalize(allvals, 'range');
                    elseif ~obj.posvisible && obj.negvisible
                        cmap = ea_colorgradient(gradientLevel, linecols(group,:), [1,1,1]);
                        voxcmap{group} = ea_colorgradient(gradientLevel, linecols(group,:), cmap(shiftedCmapEnd,:));
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
                        voxcmap{group} = [cmapLeft;cmapRight];
                        cmapind = ones(size(allvals))*gradientLevel/2;
                        cmapind(allvals<0) = round(normalize(allvals(allvals<0),'range',[1,gradientLevel/2]));
                        cmapind(allvals>0) = round(normalize(allvals(allvals>0),'range',[gradientLevel/2+1,gradientLevel]));
                        alphaind = ones(size(allvals));
                        % alphaind(allvals<0) = normalize(-1./(1+exp(-allvals(allvals<0))), 'range');
                        % alphaind(allvals>0) = normalize(1./(1+exp(-allvals(allvals>0))), 'range');
                    elseif obj.posvisible
                        cmap = ea_colorgradient(gradientLevel, [1,1,1], obj.poscolor);
                        voxcmap{group} = ea_colorgradient(gradientLevel, cmap(shiftedCmapStart,:), obj.poscolor);
                        cmapind = round(normalize(allvals,'range',[1,gradientLevel]));
                        alphaind = ones(size(allvals));
                        % alphaind = normalize(1./(1+exp(-allvals)), 'range');
                    elseif obj.negvisible
                        cmap = ea_colorgradient(gradientLevel, obj.negcolor, [1,1,1]);
                        voxcmap{group} = ea_colorgradient(gradientLevel, obj.negcolor, cmap(shiftedCmapEnd,:));
                        cmapind = round(normalize(allvals,'range',[1,gradientLevel]));
                        alphaind = ones(size(allvals));
                        % alphaind = normalize(-1./(1+exp(-allvals)), 'range');
                    end
                end
                setappdata(obj.resultfig, ['voxcmap',obj.ID], voxcmap);
                
                cmapind = mat2cell(cmapind, numel(vals{group}))';
                alphaind = mat2cell(alphaind, numel(vals{group}))';
                
                
                switch obj.vizmode
                    case 'Regions'
                        % Plot voxels if any survived
                            if obj.posvisible
                                % plot positives:
                                posvox=space;
                                posvox.img(:)=0;
                                posvox.img(usedidx{group}(vals{group}>0))=vals{group}(vals{group}>0);
                                
                                pobj.nii=posvox;
                                pobj.name='Positive';
                                pobj.niftiFilename='Positive.nii';
                                pobj.binary=0;
                                pobj.usesolidcolor=0;
                                pobj.color=obj.poscolor;
                                pobj.colormap=ea_colorgradient(length(gray), [1,1,1], obj.poscolor);
                                pobj.smooth=10;
                                pobj.hullsimplify=0.5;
                                obj.drawobject{group}{1}=ea_roi('Positive.nii',pobj);
                                
                            end
                            
                            if obj.negvisible
                                % plot negatives:
                                negvox=space;
                                negvox.img(:)=0;
                                negvox.img(usedidx{group}(vals{group}<0))=-vals{group}(vals{group}<0);
                                
                                pobj.nii=negvox;
                                pobj.name='Negative';
                                pobj.niftiFilename='Negative.nii';
                                pobj.binary=0;
                                pobj.usesolidcolor=0;
                                pobj.color=obj.negcolor;
                                pobj.colormap=ea_colorgradient(length(gray), obj.negcolor, [1,1,1]);
                                pobj.smooth=10;
                                pobj.hullsimplify=0.5;
                                obj.drawobject{group}{2}=ea_roi('Negative.nii',pobj);
                            end
                        
                    case 'Surface (Elvis)'
                        % first draw correct surface
                        switch obj.model
                            case 'Smoothed'
                                if obj.modelRH; rh=ea_readObj([ea_space,'surf_r_smoothed.obj']); end
                                if obj.modelLH; lh=ea_readObj([ea_space,'surf_l_smoothed.obj']); end
                            case 'Full'
                                if obj.modelRH; rh=ea_readObj([ea_space,'surf_r.obj']); end
                                if obj.modelLH; lh=ea_readObj([ea_space,'surf_l.obj']); end
                        end
                        cmap=ea_colorgradient(256,obj.negcolor,[1,1,1],obj.poscolor);
                        % get colors for surface:
                        bb=space.mat*[1,size(space.img,1);1,size(space.img,2);1,size(space.img,3);1,1];
                        [X,Y,Z]=meshgrid(linspace(bb(2,1),bb(2,2),size(space.img,2)),...
                            linspace(bb(1,1),bb(1,2),size(space.img,1)),...
                            linspace(bb(3,1),bb(3,2),size(space.img,3)));
                        space.img(:)=vals{group};
                        if ~obj.posvisible
                            space.img(space.img>0)=0;
                        end
                        if ~obj.negvisible
                            space.img(space.img<0)=0;
                        end
                        if obj.modelRH
                            rh_nc=round(ea_contrast(isocolors(X,Y,Z,space.img,rh.vertices))*255+1);
                            rh_nc(isnan(rh_nc))=128; % set to white for now
                            obj.drawobject{group}{1}=patch('Faces',rh.faces,'Vertices',rh.vertices,'FaceColor','interp','EdgeColor','none','FaceVertexCData',cmap(rh_nc,:),...
                                'SpecularStrength',0.35,'SpecularExponent',30,'SpecularColorReflectance',0,'AmbientStrength',0.07,'DiffuseStrength',0.45,'FaceLighting','gouraud');
                        end
                        if obj.modelLH
                            lh_nc=round(ea_contrast(isocolors(X,Y,Z,space.img,lh.vertices))*255+1);
                            lh_nc(isnan(lh_nc))=128; % set to white for now
                            obj.drawobject{group}{2}=patch('Faces',lh.faces,'Vertices',lh.vertices,'FaceColor','interp','EdgeColor','none','FaceVertexCData',cmap(lh_nc,:),...
                                'SpecularStrength',0.35,'SpecularExponent',30,'SpecularColorReflectance',0,'AmbientStrength',0.07,'DiffuseStrength',0.45,'FaceLighting','gouraud');
                        end
                    case 'Surface (Surfice)'
                        keyboard
                        
                end
                
                
                setappdata(obj.resultfig,['dt_',obj.ID],obj.drawobject); % store handle of surf to figure.
                
                
                % Set colorbar tick positions and labels
                if ~isempty(allvals)
                    if obj.posvisible && obj.negvisible
                        tick{group} = [1, length(voxcmap{group})];
                        poscbvals = sort(allvals(allvals>0));
                        negcbvals = sort(allvals(allvals<0));
                        ticklabel{group} = [negcbvals(1), poscbvals(end)];
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
            
            % store colorbar in object
            if exist('fibcmap','var') % could be no fibers present at all.
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
