function [vals,fibcell,usedidx] = ea_explorer_calcstats(obj,patientselection,OutcomePerm)

% NB: for PCA, we are going to reassign Outcome later in the function
if ~exist('OutcomePerm','var')
    Outcome=obj.responsevar;
else % used in permutation based statistics - in this case the real improvement can be substituted with permuted variables.
    Outcome=OutcomePerm;
end
if ~exist('patientselection','var') % patientselection can be supplied directly (in this case, obj.patientselection is ignored), e.g. for cross-validations.
    patientselection=obj.patientselection;
end

% fiber values can be sigmoid transformerd
if obj.SigmoidTransform 
    myvals_raw = obj.results.(ea_conn2connid(obj.connectome)).(ea_method2methodid(obj)).fibsval;
    myvals = myvals_raw;  % initialize
    for side = 1:size(myvals_raw,2)
        myvals{1,side}(:,:) = ea_SigmoidFromEfield(myvals_raw{1,side}(:,:));
    end
else
    myvals = cellfun(@full, obj.results.(ea_conn2connid(obj.connectome)).(ea_method2methodid(obj)).fibsval, 'Uni', 0);
end


if size(Outcome,2)==1 % 1 entry per patient, not per electrode
    Outcome=[Outcome,Outcome]; % both sides the same;
end

if obj.mirrorsides
    Outcome=[Outcome;Outcome];
end

if ~isempty(obj.covars)
    for Outcome=1:length(obj.covars)
        if obj.mirrorsides
            covars{Outcome} = [obj.covars{Outcome};obj.covars{Outcome}];
        else
            covars{Outcome} = obj.covars{Outcome};
        end
    end
end

if strcmp(obj.multitractmode,'Split & Color By PCA')
    % prep PCA: here, we need to get scores for all the patients. But PCA should
    % remain based on the patients selected for analysis 
    
    % variables for the selected patients 
    for Outcome=1:length(obj.subscore.vars)
        selected_subscores{Outcome} = obj.subscore.vars{Outcome}(obj.patientselection);
    end
    selected_subscores = cell2mat(selected_subscores);
    subvars=ea_nanzscore(selected_subscores);
    
    % PCA - get PC scores for selected patients
    [coeff,score,latent,tsquared,explained,mu]=pca(subvars,'rows','complete');
    % [coeff,score,latent,tsquared,explained,mu]=pca(subvars,'rows','pairwise'); %pca
    
    if isempty(score) || ea_isnan(score,'any')
        score=nan(length(obj.responsevar),obj.numpcs);
    end

    obj.subscore.pcavars=cell(obj.numpcs,1);

    for pc=1:obj.numpcs
        obj.subscore.pcavars{pc}(obj.patientselection,1)=score(:,pc); %pca variables -> pca components, location of first subscore is replaced by first pc
    end
    if ~isfield(obj.subscore,'pcacolors')
        obj.subscore.pcacolors=ea_color_wes('all'); % assign some random colors.
    end

    % now use PCA weights to get PC scores for the non-selected patients 
    patientnonsel = logical(ones(1, length(obj.allpatients)));
    patientnonsel(obj.patientselection) = 0; 

    for Outcome=1:length(obj.subscore.vars)
        nonsel_subscores{Outcome} = obj.subscore.vars{Outcome}(patientnonsel);
    end
    nonsel_subscores = cell2mat(nonsel_subscores);
    % pseudo zscore - use mean and sd of selected patients to keep same "scale"
    for ci = 1:size(selected_subscores, 2)
        datawonan = selected_subscores(:, ci);
        datawonan = datawonan(~isnan(datawonan));
        datamean(ci) = mean(datawonan);
        datasd(ci) = std(datawonan);
    end
    nonsel_subvars = ( nonsel_subscores - repmat(datamean, size(nonsel_subscores, 1), 1) ) ...
        ./ repmat(datasd, size(nonsel_subscores, 1), 1);
    
    % multiply clinical scores by weights
    for pc=1:obj.numpcs
        obj.subscore.pcavars{pc}(patientnonsel,1)= nonsel_subvars*coeff(:,pc);
    end

    % save the PCA coefficients for later 
    obj.subscore.pcacoeff = coeff; 

end

switch obj.multitractmode
    case 'Split & Color By Group'
        groups = unique(obj.M.patient.group)';
        dogroups = 1;
        dosubscores = 0;
    case 'Split & Color By Subscore'
        if ~isempty(obj.subscore.vars) %this will be empty when user
            %initializes the split by subscore button
           groups = 1:length(obj.subscore.vars);
           dosubscores = 1;
           dogroups = 0;
        else
            obj.multitractmode = 'Single Tract Analysis';
            groups = 1;
            dogroups = 0;
            dosubscores = 0;
        end
    case 'Split & Color By PCA'
        groups = 1:length(obj.subscore.pcavars);
        dosubscores = 1;
        dogroups = 0;
    otherwise
        groups=1;
        dogroups = 0;
        dosubscores = 0;
end

for group=groups
    gmyvals=myvals; %refresh myvals
    if dogroups
        groupspt=find((obj.M.patient.group)'==group);
        gpatsel=groupspt(ismember(groupspt,patientselection));
    elseif dosubscores
        gpatsel=patientselection;
    else
        gpatsel=patientselection;
    end
    if obj.mirrorsides
        gpatsel=[gpatsel,gpatsel + length(obj.allpatients)];
    end
    if dosubscores
        switch obj.multitractmode
            case 'Split & Color By Subscore'
                Outcome = obj.subscore.vars{group};
            case 'Split & Color By PCA'
                if ~exist('OutcomePerm','var')
                    Outcome = obj.subscore.pcavars{group};
                else 
                    % recompute permuted PC scores based on permuted
                    % clinical vars - keep right dimensions for Outcome but use
                    % only selected patients for zscores - patientsel
                    % applied to Outcome later
                    Outcome = nan(size(OutcomePerm)); 
                    Outcome(obj.patientselection,:) = ea_nanzscore(OutcomePerm(obj.patientselection,:))*coeff; 
                    Outcome = Outcome(:, group); 
                end    
        end
        if size(Outcome,2)==1 % 1 entry per patient, not per electrode
            Outcome=[Outcome,Outcome]; % both sides the same;
        end
        if obj.mirrorsides
            Outcome=[Outcome;Outcome];
        end
    end

    for side=1:numel(gmyvals)
        % check connthreshold
        if obj.runwhite || strcmp(obj.statmetric,'Plain Connections')
            sumgmyvals=sum(gmyvals{side}(:,gpatsel),2);
        else
            switch obj.statmetric
                case {'Two-Sample T-Tests / VTAs (Baldermann 2019) / PAM (OSS-DBS)','One-Sample Tests / VTAs / PAM (OSS-DBS)','Proportion Test (Chi-Square) / VTAs (binary vars)','Binomial Tests / VTAs (binary vars)'}
                    sumgmyvals=sum(gmyvals{side}(:,gpatsel),2);
                case {'Correlations / E-fields (Irmen 2020)','Reverse T-Tests / E-Fields (binary vars)','Odds Ratios / EF-Sigmoid (Jergas 2023)','Weighted Linear Regression / EF-Sigmoid (Dembek 2023)'}
                    if obj.SigmoidTransform == 1 && (strcmp(ea_method2methodid(obj), 'spearman_5peak') || strcmp(ea_method2methodid(obj), 'spearman_peak'))
                        % 0.5 V / mm -> 0.5 probability 
                        sumgmyvals=sum((gmyvals{side}(:,gpatsel)>obj.efieldthreshold/1000.0),2);
                    else
                        sumgmyvals=sum((gmyvals{side}(:,gpatsel)>obj.efieldthreshold),2);
                    end
            end
        end
        % remove fibers that are not connected to enough VTAs/Efields or connected
        % to too many VTAs (connthreshold slider)
        if ~obj.runwhite
            gmyvals{side}(sumgmyvals<((obj.connthreshold/100)*length(gpatsel)),gpatsel)=0;
            if ~(ismember(obj.statmetric,{'Correlations / E-fields (Irmen 2020)','Reverse T-Tests / E-Fields (binary vars)'})) % efields & reverse t-tests for binary vars cases
                % only in case of VTAs (given two-sample-t-test statistic) do we
                % need to also exclude if tract is connected to too many VTAs:
                gmyvals{side}(sumgmyvals>((1-(obj.connthreshold/100))*length(gpatsel)),gpatsel)=0;
            end
        end
        
        % init outputvars
        vals{group,side}=nan(size(gmyvals{side},1),1);
        if obj.visualization.showsignificantonly
            pvals{group,side}=vals{group,side};
        end
        if obj.runwhite || strcmp(obj.statmetric,'Plain Connections')
            vals{group,side} = sumgmyvals/length(gpatsel);
       
            if ~obj.runwhite % the white fibers will always show connection to any single vta/roi.
                vals{group,side}(sumgmyvals<((obj.connthreshold/100)*length(gpatsel)))=nan;
            else
                vals{group,side}(vals{group,side}==0)=nan;
            end

            if strcmp(obj.statmetric,'Plain Connections') && obj.visualization.showsignificantonly
                ea_error('Calculating significance does not make sense (plain connections mode)');
            end
        else
            switch obj.statmetric
                case 'Two-Sample T-Tests / VTAs (Baldermann 2019) / PAM (OSS-DBS)' % two-sample t-tests / OSS-DBS
                    % check if covariates exist:
                    if exist('covars', 'var')
                        % they do:
                        nixfib=find(any(gmyvals{side}(:,gpatsel)').*~all(gmyvals{side}(:,gpatsel)'));
                        warning('off', 'stats:classreg:regr:lmeutils:StandardLinearMixedModel:Message_PerfectFit');
                        ea_dispercent(0,['Side ',num2str(side),': Calculating T-values (taking covariates into account)']);
                        for fib=1:length(nixfib)
                            data = table(Outcome(gpatsel,side),gmyvals{side}(nixfib(fib),gpatsel)',...
                                'VariableNames',{'response','myvals'});
                            formula = 'response ~ 1 + myvals';
                            data.myvals=categorical(data.myvals);
                            for cv=1:length(covars)
                                thiscv=covars{cv}(gpatsel,:);
                                if (size(thiscv,2)==2)
                                    thiscv=thiscv(:,side);
                                end
                                data=addvars(data,thiscv,'NewVariableNames',ea_space2sub(obj.covarlabels{cv}));
                                if isequal(thiscv,logical(thiscv))
                                    formula=[formula,' + (1 + myvals | ',ea_space2sub(obj.covarlabels{cv}),')'];
                                    data.(ea_space2sub(obj.covarlabels{cv}))=categorical(data.(ea_space2sub(obj.covarlabels{cv})));
                                else
                                    formula=[formula,' + ',ea_space2sub(obj.covarlabels{cv})];
                                end
                            end
                            mdl = fitlme(data,formula);
                            vals{group,side}(nixfib(fib))=mdl.Coefficients.tStat(2);
                            if obj.visualization.showsignificantonly
                                pvals{group,side}(nixfib(fib))=mdl.Coefficients.pValue(2);
                            end
                            ea_dispercent(fib/length(nixfib));
                        end
                        ea_dispercent(1,'end');
                        fprintf('\b');
                        warning('on', 'stats:classreg:regr:lmeutils:StandardLinearMixedModel:Message_PerfectFit');

                    else
                        % no covariates exist:
                        allvals=repmat(Outcome(gpatsel,side)',size(gmyvals{side}(:,gpatsel),1),1); % improvement values (taken from Lead group file or specified in line 12).
                        fibsimpval=allvals; % Make a copy to denote improvements of connected fibers
                        fibsimpval=double(fibsimpval); % in case entered logical
                        fibsimpval(~logical(gmyvals{side}(:,gpatsel)))=nan; % Delete all unconnected values
                        nfibsimpval=allvals; % Make a copy to denote improvements of unconnected fibers
                        nfibsimpval(logical(gmyvals{side}(:,gpatsel)))=nan; % Delete all connected values
                        [~,ps,~,stats]=ttest2(fibsimpval',nfibsimpval'); % Run two-sample t-test across connected / unconnected values
                        vals{group,side}=stats.tstat';
                        
                        %% Comment this if you don't want the graphs
                        %% if you need negative fibers just use min here 
                        %[~, maxidx] = max(vals{group,side}); 
                        %imp_conn = fibsimpval(maxidx, :); 
                        %imp_nonconn = nfibsimpval(maxidx, :); 
                        %% [~,ps,~,stats]=ttest2(imp_conn', imp_nonconn'); 
                        %figure
                        %%swarmchart([ones(size(imp_conn)); ones(size(imp_conn))*2]',[imp_conn; imp_nonconn]'); 
                        %boxplot([imp_conn; imp_nonconn]');
                        %xticks(1:2)
                        %xticklabels({'Connected VTAs', 'Non-connected VTAs'}); 
                        %ylabel('Improvement')                        
                        
                        if obj.visualization.showsignificantonly
                            pvals{group,side}=ps';
                        end
                        
                        %vals{group,side}(p>0.5)=nan; % discard noisy fibers (optional or could be adapted)
                    end
                case 'One-Sample Tests / VTAs / PAM (OSS-DBS)'
                    switch obj.corrtype
                        case 'T-Tests'
                            if exist('covars', 'var')
                                ea_error('Covariates not implemented for One-Sample T-Tests.')
                            else
                                % no covariates exist:
                                allvals=repmat(Outcome(gpatsel,side)',size(gmyvals{side}(:,gpatsel),1),1); % improvement values (taken from Lead group file or specified in line 12).
                                fibsimpval=allvals; % Make a copy to denote improvements of connected fibers
                                fibsimpval=double(fibsimpval); % in case entered logical
                                fibsimpval(~logical(gmyvals{side}(:,gpatsel)))=nan; % Delete all unconnected values
                                [~,ps,~,stats]=ttest(fibsimpval'); % Run one-sample t-test across connected / unconnected values
                                vals{group,side}=stats.tstat';
                                if obj.visualization.showsignificantonly
                                    pvals{group,side}=ps';
                                end
                            end
                    
                        case 'Wicoxon Signed Rank Tests'
                            if exist('covars', 'var')
                                ea_error('Covariates not implemented for Wilcoxon Tests.')
                            else
                                ea_error('Wilcoxon Tests not yet implemented.')
                            end

                    end
                    
                case 'Correlations / E-fields (Irmen 2020)'
                    if ismember(lower(obj.corrtype),{'pearson','spearman'})
                        conventionalcorr=1;
                    else
                        conventionalcorr=0;
                    end

                    nonempty=sum(gmyvals{side}(:,gpatsel),2)>0;
                    invals=gmyvals{side}(nonempty,gpatsel)';
                    if ~isempty(invals)
                        if exist('covars', 'var') && conventionalcorr % partial corrs only implemented for Pearson & Spearman
                            usecovars=[];

                            for cv=1:length(covars)
                                thiscv=covars{cv}(gpatsel,:);
                                if (size(thiscv,2)==2)
                                    thiscv=thiscv(:,side);
                                end
                                usecovars=[usecovars,thiscv];
                            end
                            if obj.visualization.showsignificantonly
                                [outvals,outps]=partialcorr(invals,Outcome(gpatsel,side),usecovars,'rows','pairwise','type',obj.corrtype); % generate optimality values on all but left out patients
                            else % no need to calc p-val here
                                outvals=partialcorr(invals,Outcome(gpatsel,side),usecovars,'rows','pairwise','type',obj.corrtype); % generate optimality values on all but left out patients
                            end
                        else
                            if conventionalcorr
                                if obj.visualization.showsignificantonly
                                    [outvals,outps]=corr(invals,Outcome(gpatsel,side),'rows','pairwise','type',obj.corrtype); % generate optimality values on all but left out patients
                                else % no need to calc p-val here
                                    outvals=corr(invals,Outcome(gpatsel,side),'rows','pairwise','type',obj.corrtype); % generate optimality values on all but left out patients
                                end
                            else
                                if exist('covars', 'var')
                                    ea_error('Inclusion of covariates not implemented for this type of correlation');
                                end
                                switch lower(obj.corrtype)
                                    case 'bend'
                                        if obj.visualization.showsignificantonly
                                            [outvals,outps]=ea_bendcorr(invals,Outcome(gpatsel,side)); % generate optimality values on all but left out patients
                                        else % no need to calc p-val here
                                            outvals=ea_bendcorr(invals,Outcome(gpatsel,side)); % generate optimality values on all but left out patients
                                        end
                                    case 'skipped pearson'
                                        if obj.visualization.showsignificantonly
                                            ea_error('Significance not implemented for Skipped correlations');
                                        else
                                            outvals=ea_skipped_correlation(invals,Outcome(gpatsel,side),'Pearson'); % generate optimality values on all but left out patients
                                        end
                                    case 'skipped spearman'
                                        if obj.visualization.showsignificantonly
                                            ea_error('Significance not implemented for Skipped correlations');
                                        else
                                            outvals=ea_skipped_correlation(invals,Outcome(gpatsel,side),'Spearman'); % generate optimality values on all but left out patients
                                        end
                                end

                            end
                        end

                        vals{group,side}(nonempty)=outvals;
                        if exist('outps','var') % only calculated if testing for significance.
                            pvals{group,side}(nonempty)=outps;
                        end
                    end
                case 'Proportion Test (Chi-Square) / VTAs (binary vars)'

                    nonempty=sum(gmyvals{side}(:,gpatsel),2)>0;
                    invals=gmyvals{side}(nonempty,gpatsel)';
                    if ~isempty(invals)

                        ImpBinary=double((Outcome(gpatsel,side))>0); % make sure variable is actually binary
                        % restore nans
                        ImpBinary(isnan(Outcome(gpatsel,side)))=nan;
                        suminvals=sum(invals(ImpBinary == 1,:),1); % for each fiber, how many vtas cover it of patients that also had the effect (binary outcome)
                        Ninvals=sum(ImpBinary == 1,1);
                        sumoutvals=sum(invals(ImpBinary == 0,:),1); % for each fiber, how many vtas cover it of patients that did not have the effect (binary var)
                        Noutvals=sum(ImpBinary == 0,1);

                        prop=zeros(size(invals,2),1); %
                        outps=prop; %

                        for fib=1:size(invals,2)
                            [h,outps(fib), prop(fib)]  = ea_prop_test([suminvals(fib),sumoutvals(fib)],[Ninvals,Noutvals],1);
                        end

                        vals{group,side}(nonempty)=prop;

                        if exist('outps','var') % only calculated if testing for significance.
                            pvals{group,side}(nonempty)=outps;
                        end
                    end
                case 'Binomial Tests / VTAs (binary vars)'
                    nonempty=sum(gmyvals{side}(:,gpatsel),2)>0; % number of connected tracts
                    invals=gmyvals{side}(nonempty,gpatsel)';
                    if ~isempty(invals)

                        ImpBinary=double((Outcome(gpatsel,side))>0); % make sure variable is actually binary
                        % restore nans
                        ImpBinary(isnan(Outcome(gpatsel,side)))=nan;
                        sum_coverage_createeffect=sum(invals(ImpBinary == 1,:),1); % sum of coverage that creates the effect for all tracts
                        sum_coverage_noeffect=sum(invals(ImpBinary == 0,:),1); % sum of coverage that does not create the effect for all tracts

                        prob_createeffect=ea_nansum(ImpBinary)/sum(~isnan(ImpBinary));
                        sum_coverage=sum_coverage_noeffect+sum_coverage_createeffect;
                        outps = binopdf(sum_coverage_createeffect,sum_coverage_noeffect+sum_coverage_createeffect,ea_nansum(ImpBinary)/sum(~isnan(ImpBinary)));
                        thisvals = (sum_coverage_createeffect./(sum_coverage)) - prob_createeffect;
                        
                        vals{group,side}(nonempty)=thisvals;
                        if exist('outps','var') % only calculated if testing for significance.
                            pvals{group,side}(nonempty)=outps;
                        end

                    end
                case 'Reverse T-Tests / E-Fields (binary vars)'
                    nonempty=sum(gmyvals{side}(:,gpatsel),2)>0;
                    invals=gmyvals{side}(nonempty,gpatsel)';
                    if ~isempty(invals)
                        ImpBinary=double((Outcome(gpatsel,side))>0); % make sure variable is actually binary
                        % restore nans
                        ImpBinary(isnan(Outcome(gpatsel,side)))=nan;
                        upSet=invals(ImpBinary==1,:);
                        downSet=invals(ImpBinary==0,:);

                        if obj.visualization.showsignificantonly
                            [~,ps,~,stats]=ttest2(upSet,downSet); % Run two-sample t-test across connected / unconnected values
                            outvals=stats.tstat';
                            outps=ps;
                        else % no need to calc p-val here
                            [~,~,~,stats]=ttest2(upSet,downSet); % Run two-sample t-test across connected / unconnected values
                            outvals=stats.tstat';
                        end

                        vals{group,side}(nonempty)=outvals;
                        if exist('outps','var') % only calculated if testing for significance.
                            pvals{group,side}(nonempty)=outps;
                        end
                    end
                case 'Weighted Linear Regression / EF-Sigmoid (Dembek 2023)'
                    nonempty=sum(gmyvals{side}(:,gpatsel),2)>0;
                    invals=gmyvals{side}(nonempty,gpatsel)';
    
                    %weighted_linear_regression_test = 'two-sample';
                    %if strcmp('two-sample',weighted_linear_regression_test)
                    if strcmp('2-Sample',obj.twoSampleWeighted)
                        [outvals,outps] = ea_discfibers_TwoSample_weightedLinearRegression(invals,Outcome(gpatsel,side),'t-value');
                    else
                        [outvals,outps] = ea_discfibers_weightedLinearRegression(invals,Outcome(gpatsel,side),'mean','t-value'); % generate optimality values on all but left out patients
                    end
                    vals{group,side}(nonempty)=outvals;
                    if exist('outps','var') % only calculated if testing for significance.
                        pvals{group,side}(nonempty)=outps;
                    end
                case 'Odds Ratios / EF-Sigmoid (Jergas 2023)'

                    nonempty=sum(gmyvals{side}(:,gpatsel),2)>0; % number of connected tracts
                    invals=gmyvals{side}(nonempty,gpatsel)';

                    if ~isempty(invals)
                    
                        Impr = Outcome(gpatsel,side);

                        if any(isnan(Outcome(gpatsel,side)))
                            ea_warndlg("NaNs in scores detected, removing...")
                            %return
                            % remove the whole row
                            nan_scores = isnan(Outcome(gpatsel,side));
                            
                            Impr(nan_scores,:) = [];
                            invals(nan_scores,:) = [];
                        end

                        ImpBinary=logical((Impr)>0); % make sure variable is actually binary
                        % restore nans

                        [outvals,CI95_up,CI95_low,outps] = ea_discfibers_odds_ratios(invals,ImpBinary);
             
                        vals{group,side}(nonempty)=outvals;
                        if exist('outps','var') % only calculated if testing for significance.
                            pvals{group,side}(nonempty)=outps;
                        end
                    end
                    
            end
        end
    end
end

% close group loop to test for significance across all tests run:
if ~obj.runwhite
    if obj.visualization.showsignificantonly
        vals=ea_corrsignan(vals,pvals,obj);
    end
end

% reopen group loop for thresholding etc:
for group=groups
    for side=1:numel(gmyvals)
        fibcell{group,side}=obj.results.(ea_conn2connid(obj.connectome)).fibcell{side}(~isnan(vals{group,side}));
        % Remove vals and fibers outside the thresholding range
        obj.stats.pos.available(side)=sum(cat(1,vals{:,side})>0); % only collected for first group (positives)
        obj.stats.neg.available(side)=sum(cat(1,vals{:,side})<0);
        if dosubscores || dogroups
            if ~obj.subscore.special_case
                obj.subscore.vis.pos_available(group,side)=sum(cat(1,vals{group,side})>0); % collected for every group
                obj.subscore.vis.neg_available(group,side)=sum(cat(1,vals{group,side})<0);
            end
        end
        usedidx{group,side}=find(~isnan(vals{group,side}));
        vals{group,side}=vals{group,side}(usedidx{group,side}); % final weights for surviving fibers
        if exist('pvals','var')
            pvals{group,side}=pvals{group,side}(usedidx{group,side}); % final weights for surviving fibers
        end


        switch obj.threshstrategy
            case 'Fixed Amount' % here we want to create threshs for each side separately.
                posvals = sort(vals{group,side}(vals{group,side}>0),'descend');
                negvals = sort(vals{group,side}(vals{group,side}<0),'ascend');
            otherwise % in other cases, we want to apply the same thresh to both sides.
                allvals = vertcat(vals{group,:});
                posvals = sort(allvals(allvals>0),'descend');
                negvals = sort(allvals(allvals<0),'ascend');
        end
        % positive thresholds
        if dosubscores || dogroups
            if obj.subscore.special_case
                if ~obj.visualization.posvisible || ~obj.visualization.showposamount(side) || isempty(posvals)
                    posthresh = inf;
                else
                    posthresh = ea_fibValThresh(obj.threshstrategy, posvals, obj.visualization.showposamount(side));
                end
            else
                if ~obj.subscore.posvisible(group) || ~obj.subscore.vis.showposamount(group,side) || isempty(posvals)
                    posthresh = inf;
                else
                    posthresh = ea_fibValThresh(obj.threshstrategy, posvals, obj.subscore.vis.showposamount(group,side));
                end
            end
        else
            if ~obj.visualization.posvisible || ~obj.visualization.showposamount(side) || isempty(posvals)
                posthresh = inf;
            else
                posthresh = ea_fibValThresh(obj.threshstrategy, posvals, obj.visualization.showposamount(side));
            end
        end

        % negative thresholds
        if dosubscores || dogroups
            if obj.subscore.special_case
                if ~obj.visualization.negvisible || ~obj.visualization.shownegamount(side) || isempty(negvals)
                    negthresh = -inf;
                else
                    negthresh = ea_fibValThresh(obj.threshstrategy, negvals, obj.visualization.shownegamount(side));
                end
            else
                if ~obj.subscore.negvisible(group) || ~obj.subscore.vis.shownegamount(group,side) || isempty(negvals)
                    negthresh = -inf;
                else
                    negthresh = ea_fibValThresh(obj.threshstrategy, negvals, obj.subscore.vis.shownegamount(group,side));
                end
            end
        else
            if ~obj.visualization.negvisible || ~obj.visualization.shownegamount(side) || isempty(negvals)
                negthresh = -inf;
            else
                negthresh = ea_fibValThresh(obj.threshstrategy, negvals, obj.visualization.shownegamount(side));
            end
        end
        if ~obj.runwhite
            % Remove vals and fibers outside the thresholding range (set by
            % sliders)
            remove = logical(logical(vals{group,side}<posthresh) .* logical(vals{group,side}>negthresh));
            vals{group,side}(remove)=[];
            fibcell{group,side}(remove)=[];
            usedidx{group,side}(remove)=[];
        end
    end
end



function vals=ea_corrsignan(vals,ps,obj)
allvals=cat(1,vals{:});
allps=cat(1,ps{:});

nnanidx=~isnan(allvals);
numtests=sum(nnanidx);

switch lower(obj.multcompstrategy)
    case 'fdr'
        pnnan=allps(nnanidx); % pvalues from non-nan value entries
        [psort,idx]=sort(pnnan);
        pranks=zeros(length(psort),1);
        for rank=1:length(pranks)
            pranks(idx(rank))=rank;
        end
        pnnan=pnnan.*numtests;
        pnnan=pnnan./pranks;
        allps(nnanidx)=pnnan;
    case 'bonferroni'
        allps(nnanidx)=allps(nnanidx).*numtests;
end

allps(~nnanidx)=1; % set p values from values with nan to 1
allvals(allps>obj.alphalevel)=nan; % delete everything nonsignificant.

% feed back into cell format:
cnt=1;
for cellentry=1:numel(vals)
    vals{cellentry}(:)=allvals(cnt:cnt+length(vals{cellentry})-1);
    ps{cellentry}(:)=allps(cnt:cnt+length(vals{cellentry})-1);
    cnt=cnt+length(vals{cellentry});
end


function fibValThreshold = ea_fibValThresh(threshstrategy, vals, threshold)
switch threshstrategy
    case 'Percentage Relative to Peak'
        range = vals(1) - vals(end);
        fibValThreshold = vals(1) - threshold/100 * range;
        if range == 0
            if vals(1) > 0
                fibValThreshold = fibValThreshold - eps*10;
            else
                fibValThreshold = fibValThreshold + eps*10;
            end
        end
    case 'Percentage Relative to Amount'
        index = round((threshold/100)*length(vals));
        if index <=0
            fibValThreshold = vals(1);
        else
            fibValThreshold = vals(index);
        end
    case 'Fixed Amount'
        if length(vals)>round(threshold)
            fibValThreshold=vals(round(threshold));
        else
            fibValThreshold=vals(end);
        end
    case 'Histogram (CDF)'
        if vals(1) > 0
            [fx, x] = ecdf(vals);
            fibValThreshold = x(find(fx>=(1-threshold), 1));
        else
            [fx, x] = ecdf(-vals);
            fibValThreshold = -x(find(fx>=(1-threshold), 1));
        end
    case 'Fixed Fiber Value'
        fibValThreshold = threshold;
end

function result = ea_isnan(input_array,flag)
if size(input_array,2) > 1
    if strcmp(flag,'any')
        op = any(isnan(input_array));
    elseif strcmp(flag,'all')
        op = all(isnan(input_array));
    end
    if length(find(op)) > 1
        result = 1;
    else
        result = 0;

    end
else
    if strcmp(flag,'any')
        result = any(isnan(input_array));
    elseif strcmp(flag,'all')
        result = all(isnan(input_array));
    end
end
