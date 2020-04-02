classdef ea_disctract < handle
    % Discriminative fiber class to handle visualizations of discriminative fibers in lead dbs resultfig / 3D Matlab figures
    % A. Horn

    properties (SetObservable)
        M % content of lead group project
        tractset % parent object of ea_disctractset
        resultfig % figure handle to plot results
        ID % name / ID of discriminative fibers object
        posvisible % pos tract visible
        negvisible % neg tract visible
        showposamount % two entries for right and left
        shownegamount % two entries for right and left
        connthreshold
        statmetric % entry from discfiber settings as initially specified in prefs.machine.lg (which are separately stored for each analysis/object).
        poscolor % positive main color
        negcolor % negative main color
        % the  main content variables:

        results
        % includes subfields by results.connectomename.ttests /
        % results.connectomename.efields with
        % fibcell % cell of all fibers connected
        % fibsval % connection weight value for each fiber to each VTA
        % fibweights % usually T- or R-values associated with each tract
        %
        patientselection % selected patients to include. Note that connected fibers are always sampled from all (& mirrored) VTAs of the lead group file
        allpatients % list of all patients (as from M.patient.list)
        mirrorsides % flag to mirror VTAs / Efields to contralateral sides using ea_flip_lr_nonlinear()
        responsevar % response variable
        responsevarlabel % label of response variable
        covars % covariates
        covarlabels % covariate labels
        leadgroup % redundancy protocol only, path to original lead group project
        connectome % redundancy protocol only, name of underlying connectome
        colorbar % colorbar information
    end

    properties (Access = private)
        switchedFromSpace=3 % if switching space, this will protocol where from
    end

    methods
        function obj=ea_disctract() % class constructor

        end

        function calculate(obj)
            % check that this has not been calculated before:
            if ~isempty(obj.results) % something has been calculated
                if isfield(obj.results,conn2connid(obj.connectome))

                    if isfield(obj.results.(conn2connid(obj.connectome)),method2methodid(obj.statmetric)) % this combination was already calculated.
                        return
                    end
                end
            end

            %%% DEBUGGING
            %d=load('/PA/Neuro/_projects/lead/lead_demo/testresult.mat');
            %obj.results=d.results;
            %return
            %%%

            efieldthresh=obj.connthreshold;

            options.native = 0;

            allroilist=cell(length(obj.allpatients)*2,2);
            switch obj.statmetric
                case 1 % use paired T-Tests and binary VTA
                    suffix='';
                case 2 % use Spearman Rs and E-Fields
                    suffix='_efield';
                    prefs=ea_prefs;
                    if strcmp(prefs.lcm.vatseed,'efield_gauss')
                        suffix='_efield_gauss';
                    end
            end
            cnt=1;
            for sub=1:length(obj.allpatients) % all patients - for connected fibers selection ? and always flip
                allroilist(cnt,:)={[obj.allpatients{sub},filesep,'stimulations',filesep,ea_nt(options),'gs_',obj.M.guid,filesep,'vat',suffix,'_right.nii'],[obj.allpatients{sub},filesep,'stimulations',filesep,ea_nt(options),'gs_',obj.M.guid,filesep,'vat',suffix,'_left.nii']};
                cnt=cnt+1;
            end
            for sub=1:length(obj.allpatients) % all patients - for connected fibers selection ? and always flip
                ea_genflippedjointnii([obj.allpatients{sub},filesep,'stimulations',filesep,ea_nt(options),'gs_',obj.M.guid,filesep,'vat',suffix,'_right.nii'],[obj.allpatients{sub},filesep,'stimulations',filesep,ea_nt(options),'gs_',obj.M.guid,filesep,'vat',suffix,'_left.nii']);
                allroilist(cnt,:)={[obj.allpatients{sub},filesep,'stimulations',filesep,ea_nt(options),'gs_',obj.M.guid,filesep,'fl_','vat',suffix,'_left.nii'],[obj.allpatients{sub},filesep,'stimulations',filesep,ea_nt(options),'gs_',obj.M.guid,filesep,'fl_','vat',suffix,'_right.nii']};
                cnt=cnt+1;
            end
            cfile=[ea_getconnectomebase('dMRI'),obj.connectome,filesep,'data.mat'];
            mirroredpatselection=[obj.patientselection,obj.patientselection+length(obj.allpatients)];
            switch obj.statmetric
                case 1 % ttests
                    [fibcell,fibsval,XYZmm,nii]=ea_discfibers_heatfibertracts(cfile,{allroilist},mirroredpatselection,{obj.variables(:,1)},obj.connthreshold/100);
                case 2 % spearmans R
                    [fibcell,fibsval,XYZmm,nii,valsmm]=ea_discfibers_heatfibertracts_corr(cfile,{allroilist},mirroredpatselection,{obj.variables(:,1)},efieldthresh);
                    obj.results.(conn2connid(obj.connectome)).(method2methodid(obj.statmetric)).valsmm=valsmm;
            end

            % Main output of results - this is all we will ever need if statmetric
            % and connectome wont change
            obj.results.(conn2connid(obj.connectome)).(method2methodid(obj.statmetric)).fibsval=fibsval;
            obj.results.(conn2connid(obj.connectome)).(method2methodid(obj.statmetric)).fibcell=fibcell;
            obj.results.(conn2connid(obj.connectome)).(method2methodid(obj.statmetric)).XYZmm=XYZmm;
            obj.results.(conn2connid(obj.connectome)).(method2methodid(obj.statmetric)).nii=nii;

        end

        function [vals]=calcstats(fibsval,I,obj,patsel)

            % quickly recalc stats:
            if ~exist('patsel','var') % patsel can be supplied directly (in this case, obj.patientselection is ignored), e.g. for cross-validations.
                patsel=obj.patientselection;
            end
            if obj.mirrorsides
                patsel=[patsel;patsel+length(obj.allpatients)];
            end

            if size(I,2)==1 % 1 entry per patient, not per electrode
                I=[I,I]; % both sides the same;
            end

            for side=1:2
                switch obj.statmetric
                    case 1 % t-tests
                        allvals=repmat(I(patsel,side),1,size(fibsval{side}(:,patsel),1)); % improvement values (taken from Lead group file or specified in line 12).
                        fibsimpval=allvals; % Make a copy to denote improvements of connected fibers
                        fibsimpval(~logical(fibsval{side}(:,patsel)))=nan; % Delete all unconnected values
                        nfibsimpval=allvals; % Make a copy to denote improvements of unconnected fibers
                        nfibsimpval(logical(fibsval{side}(:,patsel)))=nan; % Delete all connected values
                        [~,p,~,stats]=ttest2(fibsimpval,nfibsimpval); % Run two-sample t-test across connected / unconnected values
                        vals{side}=stats.tstat;
                        vals{side}(p>0.5)=nan; % discard noisy fibers (optional or could be adapted)
                    case 2 % spearmans correlations
                        nanfibsval=fibsval{side}(:,patsel);
                        nanfibsval(nanfibsval==0)=nan; % only used in spearmans correlations
                        vals{side}=corr(nanfibsval',I(patsel,side),'rows','pairwise','type','Spearman'); % generate optimality values on all but left out patients
                end
            end

        end

        function draw(obj)
            I=obj.responsevar; % need to correct for bilateral vars

            [vals]=calcstats(obj.results.(conn2connid(obj.connectome)).(method2methodid(obj.statmetric)).fibsval,I,obj);


            % Contruct default blue to red colormap
            colormap(gray);
            fibcmap = ea_colorgradient(1024, [0,0,1], [1,1,1], [1,0,0]);
            setappdata(obj.resultfig, ['fibcmap',obj.ID], fibcmap);

            % Set alphas of fibers with light color to 0
            colorbarThreshold = 0.60; % Percentage of the pos/neg color to be kept
            negUpperBound=ceil(size(fibcmap,1)/2*colorbarThreshold);
            poslowerBound=floor((size(fibcmap,1)-size(fibcmap,1)/2*colorbarThreshold));

            for side=1:2
                fibcell{side}=obj.results.(conn2connid(obj.connectome)).(method2methodid(obj.statmetric)).fibcell(~isnan(vals{side}));

                vals{side}=vals{side}(~isnan(vals{side}))'; % final weights for surviving fibers

                posvals{side}=sort(vals{side}(vals{side}>0),'descend');
                negvals{side}=sort(vals{side}(vals{side}<0),'ascend');

                posthresh{side}=posvals{side}(round((obj.showposamount(side)/100)*length(posvals{side})));
                negthresh{side}=negvals{side}(round((obj.shownegamount(side)/100)*length(negvals{side})));

                % Remove vals and fibers outside the thresholding range
                remove=logical(logical(vals{side}<posthresh{side}) .* logical(vals{side}>negthresh{side}));
                vals{side}(remove)=[];
                fibcell{side}(remove)=[];

                tvalsRescale{side} = vals{side};
                tvalsRescale{side}(vals{side}>0) = ea_rescale(vals{side}(vals{side}>0), [0 1]);
                tvalsRescale{side}(vals{side}<0) = ea_rescale(vals{side}(vals{side}<0), [-1 0]);

                fibcolorInd{side}=tvalsRescale{side}*(size(fibcmap,1)/2-0.5);
                fibcolorInd{side}=fibcolorInd{side}+(size(fibcmap,1)/2+0.5);

                alphas{side}=zeros(size(fibcolorInd{side},1),1);
                if obj.posvisible && ~obj.negvisible
                	alphas{side}(round(fibcolorInd{side})>=poslowerBound) = 1;
                elseif ~obj.posvisible && obj.negvisible
                	alphas{side}(round(fibcolorInd{side})<=negUpperBound) = 1;
                elseif obj.posvisible && obj.negvisible
                    alphas{side}(round(fibcolorInd{side})>=poslowerBound) = 1;
                    alphas{side}(round(fibcolorInd{side})<=negUpperBound) = 1;
                end

                alphas{side}(round(fibcolorInd{side})>=poslowerBound) = 1;
                fibalpha=mat2cell(alphas{side},ones(size(fibcolorInd{side},1),1));

                % Plot fibers
                h=streamtube(fibcell{side},0.2);
                nones=repmat({'none'},size(fibcolorInd{side}));
                [h.EdgeColor]=nones{:};

                % Calulate fiber colors
                colors=fibcmap(round(fibcolorInd{side}),:);
                fibcolor=mat2cell(colors,ones(size(fibcolorInd{side})));

                % Set fiber colors and alphas
                [h.FaceColor]=fibcolor{:};
                [h.FaceAlpha]=fibalpha{:};
            end

            % Set colorbar tick positions and labels
            cbvals = [vals{1}(logical(alphas{1}));vals{2}(logical(alphas{2}))];
            % cbvals=tvalsRescale{side}(logical(alphas));
            if obj.posvisible && ~obj.negvisible
                cbmap = fibcmap(ceil(length(fibcmap)/2+0.5):end,:);
                tick = [poslowerBound, length(fibcmap)] - floor(length(fibcmap)/2) ;
                poscbvals = sort(cbvals(cbvals>0));
                ticklabel = [poscbvals(1), poscbvals(end)];
                ticklabel = arrayfun(@(x) num2str(x,'%.2f'), ticklabel, 'Uni', 0);
            elseif ~obj.posvisible && obj.negvisible
                cbmap = fibcmap(1:floor(length(fibcmap)/2-0.5),:);
                tick = [1, negUpperBound];
                negcbvals = sort(cbvals(cbvals<0));
                ticklabel = [negcbvals(1), negcbvals(end)];
                ticklabel = arrayfun(@(x) num2str(x,'%.2f'), ticklabel, 'Uni', 0);
            elseif obj.posvisible && obj.negvisible
                cbmap = fibcmap;
                tick = [1, negUpperBound, poslowerBound, length(fibcmap)];
                poscbvals = sort(cbvals(cbvals>0));
                negcbvals = sort(cbvals(cbvals<0));
                ticklabel = [min(cbvals), negcbvals(end), poscbvals(1), max(cbvals)];
                ticklabel = arrayfun(@(x) num2str(x,'%.2f'), ticklabel, 'Uni', 0);
            end

            % colorbar
            obj.colorbar.cmap = cbmap;
            obj.colorbar.tick = tick;
            obj.colorbar.ticklabel = ticklabel;
        end
    end

    methods (Static)
        function changeevent(~,event)
            update_trajectory(event.AffectedObject,event.Source.Name);
        end
    end
end



function id=method2methodid(method)
switch method
    case 1
        id='ttests';
    case 2
        id='spearman';
end
end

function conname=conn2connid(conname)

conname=strrep(conname,' ','');
conname=strrep(conname,'_','');
conname=strrep(conname,'(','');
conname=strrep(conname,')','');
conname=strrep(conname,'-','');
end
