function [vals,fibcell,usedidx] = ea_discfibers_calcstats(obj,patsel,Iperm)

if ~exist('Iperm','var')
    I=obj.responsevar;
else % used in permutation based statistics - in this case the real improvement can be substituted with permuted variables.
    I=Iperm;
end

fibsval = full(obj.results.(ea_conn2connid(obj.connectome)).(ea_method2methodid(obj)).fibsval);

% quickly recalc stats:
if ~exist('patsel','var') % patsel can be supplied directly (in this case, obj.patientselection is ignored), e.g. for cross-validations.
    patsel=obj.patientselection;
end

if size(I,2)==1 % 1 entry per patient, not per electrode
    I=[I,I]; % both sides the same;
end

if obj.mirrorsides
    I=[I;I];
end

if ~isempty(obj.covars)
    for i=1:length(obj.covars)
        if obj.mirrorsides
            covars{i} = [obj.covars{i};obj.covars{i}];
        else
            covars{i} = obj.covars{i};
        end
    end
end

if strcmp(obj.multitractmode,'Split & Color By PCA')
    % prep PCA:
    subvars=ea_nanzscore(cell2mat(obj.subscore.vars')); %standardize subscore stage
    try
        [coeff,score,latent,tsquared,explained,mu]=pca(subvars,'rows','pairwise'); %pca
    catch %can fail due to presence of NaN/INF TODO
        % pca failed, likely not enough variables selected.
        score=nan(length(obj.responsevar),obj.numpcs);
    end

    if isempty(score) || isnan(score)
        score=nan(length(obj.responsevar),obj.numpcs);
    end

    %if score is empty
    obj.subscore.pcavars=cell(obj.numpcs,1);
    %end

    for pc=1:obj.numpcs
        obj.subscore.pcavars{pc}=score(:,pc); %pca variables -> pca components, location of first subscore is replaced by first pc
    end
    if ~isfield(obj.subscore,'pcacolors')
        obj.subscore.pcacolors=ea_color_wes('all'); % assign some random colors.
    end
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
    gfibsval=fibsval; %refresh fibsval
    if dogroups
        groupspt=find(obj.M.patient.group==group);
        gpatsel=groupspt(ismember(groupspt,patsel));
    elseif dosubscores
        gpatsel=patsel;
    else
        gpatsel=patsel;
    end
    if obj.mirrorsides
        gpatsel=[gpatsel,gpatsel+length(obj.allpatients)];
    end
    if dosubscores
        switch obj.multitractmode
            case 'Split & Color By Subscore'
                I = obj.subscore.vars{group};
            case 'Split & Color By PCA'
                I = obj.subscore.pcavars{group};
        end
        if size(I,2)==1 % 1 entry per patient, not per electrode
            I=[I,I]; % both sides the same;
        end
        if obj.mirrorsides
            I=[I;I];
        end
    end

    for side=1:numel(gfibsval)
        % check connthreshold
        switch obj.statmetric
            case {1,3,4}
                sumgfibsval=sum(gfibsval{side}(:,gpatsel),2);
            case {2,5}
                sumgfibsval=sum((gfibsval{side}(:,gpatsel)>obj.efieldthreshold),2);
            case 6
                sumgfibsval=sum(gfibsval{side}(:,gpatsel),2);
        end
        % remove fibers that are not connected to enough VTAs/Efields or connected
        % to too many VTAs (connthreshold slider)
        gfibsval{side}(sumgfibsval<((obj.connthreshold/100)*length(gpatsel)),gpatsel)=0;
        % only in case of VTAs (given two-sample-t-test statistic) do we
        % need to also exclude if tract is connected to too many VTAs:
        if ismember(obj.statmetric,[1,4])
            gfibsval{side}(sumgfibsval>((1-(obj.connthreshold/100))*length(gpatsel)),gpatsel)=0;
        end
        % init outputvars
        vals{group,side}=nan(size(gfibsval{side},1),1);
        if obj.showsignificantonly
            pvals{group,side}=vals{group,side};
        end
        switch obj.statmetric
            case {1,3} % t-tests
                % check if covariates exist:
                if exist('covars', 'var')
                    % they do:
                    nixfib=find(any(gfibsval{side}(:,gpatsel)').*~all(gfibsval{side}(:,gpatsel)'));
                    warning('off', 'stats:classreg:regr:lmeutils:StandardLinearMixedModel:Message_PerfectFit');
                    ea_dispercent(0,['Side ',num2str(side),': Calculating T-values (taking covariates into account)']);
                    for fib=1:length(nixfib)
                        data = table(I(gpatsel,side),gfibsval{side}(nixfib(fib),gpatsel)',...
                            'VariableNames',{'response','fibsval'});
                        formula = 'response ~ 1 + fibsval';
                        data.fibsval=categorical(data.fibsval);
                        for cv=1:length(covars)
                            thiscv=covars{cv}(gpatsel,:);
                            if (size(thiscv,2)==2)
                                thiscv=thiscv(:,side);
                            end
                            data=addvars(data,thiscv,'NewVariableNames',ea_space2sub(obj.covarlabels{cv}));
                            if isequal(thiscv,logical(thiscv))
                                formula=[formula,' + (1 + fibsval | ',ea_space2sub(obj.covarlabels{cv}),')'];
                                data.(ea_space2sub(obj.covarlabels{cv}))=categorical(data.(ea_space2sub(obj.covarlabels{cv})));
                            else
                                formula=[formula,' + ',ea_space2sub(obj.covarlabels{cv})];
                            end
                        end
                        mdl = fitlme(data,formula);
                        vals{group,side}(nixfib(fib))=mdl.Coefficients.tStat(2);
                        if obj.showsignificantonly
                            pvals{group,side}(nixfib(fib))=mdl.Coefficients.pValue(2);
                        end
                        ea_dispercent(fib/length(nixfib));
                    end
                    ea_dispercent(1,'end');
                    fprintf('\b');
                    warning('on', 'stats:classreg:regr:lmeutils:StandardLinearMixedModel:Message_PerfectFit');

                else
                    % no covariates exist:
                    allvals=repmat(I(gpatsel,side)',size(gfibsval{side}(:,gpatsel),1),1); % improvement values (taken from Lead group file or specified in line 12).
                    fibsimpval=allvals; % Make a copy to denote improvements of connected fibers
                    fibsimpval(~logical(gfibsval{side}(:,gpatsel)))=nan; % Delete all unconnected values
                    nfibsimpval=allvals; % Make a copy to denote improvements of unconnected fibers
                    nfibsimpval(logical(gfibsval{side}(:,gpatsel)))=nan; % Delete all connected values
                    [~,ps,~,stats]=ttest2(fibsimpval',nfibsimpval'); % Run two-sample t-test across connected / unconnected values
                    vals{group,side}=stats.tstat';
                    if obj.showsignificantonly
                        pvals{group,side}=ps';
                    end

                    %vals{group,side}(p>0.5)=nan; % discard noisy fibers (optional or could be adapted)
                end
            case 2 % correlations
                if ismember(lower(obj.corrtype),{'pearson','spearman'})
                    conventionalcorr=1;
                else
                    conventionalcorr=0;
                end

                nonempty=full(sum(gfibsval{side}(:,gpatsel),2))>0;
                invals=gfibsval{side}(nonempty,gpatsel)';
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
                        if obj.showsignificantonly
                            [outvals,outps]=partialcorr(invals,I(gpatsel,side),usecovars,'rows','pairwise','type',obj.corrtype); % generate optimality values on all but left out patients
                        else % no need to calc p-val here
                            outvals=partialcorr(invals,I(gpatsel,side),usecovars,'rows','pairwise','type',obj.corrtype); % generate optimality values on all but left out patients
                        end
                    else
                        if conventionalcorr
                            if obj.showsignificantonly
                                [outvals,outps]=corr(invals,I(gpatsel,side),'rows','pairwise','type',obj.corrtype); % generate optimality values on all but left out patients
                            else % no need to calc p-val here
                                outvals=corr(invals,I(gpatsel,side),'rows','pairwise','type',obj.corrtype); % generate optimality values on all but left out patients
                            end
                        else
                            if exist('covars', 'var')
                                ea_error('Inclusion of covariates not implemented for this type of correlation');
                            end
                            switch lower(obj.corrtype)
                                case 'bend'
                                    if obj.showsignificantonly
                                        [outvals,outps]=ea_bendcorr(invals,I(gpatsel,side)); % generate optimality values on all but left out patients
                                    else % no need to calc p-val here
                                        outvals=ea_bendcorr(invals,I(gpatsel,side)); % generate optimality values on all but left out patients
                                    end
                                case 'skipped pearson'
                                    if obj.showsignificantonly
                                        ea_error('Significance not implemented for Skipped correlations');
                                    else
                                        outvals=ea_skipped_correlation(full(invals),I(gpatsel,side),'Pearson'); % generate optimality values on all but left out patients
                                    end
                                case 'skipped spearman'
                                    if obj.showsignificantonly
                                        ea_error('Significance not implemented for Skipped correlations');
                                    else
                                        outvals=ea_skipped_correlation(full(invals),I(gpatsel,side),'Spearman'); % generate optimality values on all but left out patients
                                    end
                            end

                        end
                    end

                    vals{group,side}(nonempty)=outvals;
                    if exist('outps','var') % only calculated if testing for significance.
                        pvals{group,side}(nonempty)=outps;
                    end
                end
            case 4 % Proportion Test (Chi-Square) / VTAs (binary vars)

                nonempty=full(sum(gfibsval{side}(:,gpatsel),2))>0;
                invals=gfibsval{side}(nonempty,gpatsel)';
                if ~isempty(invals)

                    ImpBinary=double((I(gpatsel,side))>0); % make sure variable is actually binary
                    % restore nans
                    ImpBinary(isnan(I(gpatsel,side)))=nan;
                    suminvals=full(sum(invals(logical(ImpBinary),:),1));
                    Ninvals=sum(logical(ImpBinary),1);
                    sumoutvals=full(sum(invals(~logical(ImpBinary),:),1));
                    Noutvals=sum(~logical(ImpBinary),1);

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

            case 5 % Reverse t-tests with efields for binary variables
                nonempty=full(sum(gfibsval{side}(:,gpatsel),2))>0;
                invals=gfibsval{side}(nonempty,gpatsel)';
                if ~isempty(invals)
                    ImpBinary=double((I(gpatsel,side))>0); % make sure variable is actually binary
                    % restore nans
                    ImpBinary(isnan(I(gpatsel,side)))=nan;
                    upSet=invals(ImpBinary==1,:);
                    downSet=invals(ImpBinary==0,:);

                    if obj.showsignificantonly
                        [~,ps,~,stats]=ttest2(full(upSet),full(downSet)); % Run two-sample t-test across connected / unconnected values
                        outvals=stats.tstat';
                        outps=ps;
                    else % no need to calc p-val here
                        [~,~,~,stats]=ttest2(full(upSet),full(downSet)); % Run two-sample t-test across connected / unconnected values
                        outvals=stats.tstat';
                    end

                    vals{group,side}(nonempty)=outvals;
                    if exist('outps','var') % only calculated if testing for significance.
                        pvals{group,side}(nonempty)=outps;
                    end
                end

            case 6 % Plain Connection
                vals{group,side} = sumgfibsval/length(gpatsel);
                vals{group,side}(sumgfibsval<((obj.connthreshold/100)*length(gpatsel)))=nan;
                if obj.showsignificantonly
                    ea_error('Calculating significance does not make sense (plain connections mode)');
                end
        end
    end
end

% close group loop to test for significance across all tests run:

if obj.showsignificantonly
    vals=ea_corrsignan(vals,pvals,obj);
end

% reopen group loop for thresholding etc:

for group=groups
    for side=1:numel(gfibsval)
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

        allvals = vertcat(vals{group,:});
        posvals = sort(allvals(allvals>0),'descend');
        negvals = sort(allvals(allvals<0),'ascend');
        
        % positive thresholds
        if dosubscores || dogroups
            if obj.subscore.special_case
                if ~obj.posvisible || ~obj.showposamount(side) || isempty(posvals)
                    posthresh = inf;
                else
                    posthresh = ea_fibValThresh(obj.threshstrategy, posvals, obj.showposamount(side));
                end
            else
                if ~obj.subscore.posvisible(group) || ~obj.subscore.vis.showposamount(group,side) || isempty(posvals)
                    posthresh = inf;
                else
                    posthresh = ea_fibValThresh(obj.threshstrategy, posvals, obj.subscore.vis.showposamount(group,side));
                end
            end
        else
            if ~obj.posvisible || ~obj.showposamount(side) || isempty(posvals)
                posthresh = inf;
            else
                posthresh = ea_fibValThresh(obj.threshstrategy, posvals, obj.showposamount(side));
            end
        end

        % negative thresholds
        if dosubscores || dogroups
            if obj.subscore.special_case
                if ~obj.negvisible || ~obj.shownegamount(side) || isempty(negvals)
                    negthresh = -inf;
                else
                    negthresh = ea_fibValThresh(obj.threshstrategy, negvals, obj.shownegamount(side));
                end
            else
                if ~obj.subscore.negvisible(group) || ~obj.subscore.vis.shownegamount(group,side) || isempty(negvals)
                    negthresh = -inf;
                else
                    negthresh = ea_fibValThresh(obj.threshstrategy, negvals, obj.subscore.vis.shownegamount(group,side));
                end
            end
        else
            if ~obj.negvisible || ~obj.shownegamount(side) || isempty(negvals)
                negthresh = -inf;
            else
                negthresh = ea_fibValThresh(obj.threshstrategy, negvals, obj.shownegamount(side));
            end
        end

        % Remove vals and fibers outside the thresholding range (set by
        % sliders)
        remove = logical(logical(vals{group,side}<posthresh) .* logical(vals{group,side}>negthresh));
        vals{group,side}(remove)=[];
        fibcell{group,side}(remove)=[];
        usedidx{group,side}(remove)=[];
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
        fibValThreshold = vals(round((threshold/100)*length(vals)));
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
