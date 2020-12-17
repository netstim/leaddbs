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

if obj.splitbygroup
    groups=unique(obj.M.patient.group)';
    dogroups=1;
else
    groups=1; % all one, group ignored
    dogroups=0;
end

for group=groups
    gfibsval=fibsval; %refresh fibsval
    if dogroups
        groupspt=find(obj.M.patient.group==group);
        gpatsel=groupspt(ismember(groupspt,patsel));
    else
        gpatsel=patsel;
    end
    if obj.mirrorsides
        gpatsel=[gpatsel,gpatsel+length(obj.allpatients)];
    end
    
    for side=1:numel(gfibsval)
        % check connthreshold
        switch obj.statmetric
            case 1
                sumgfibsval=sum(gfibsval{side}(:,gpatsel),2);
            case 2
                sumgfibsval=sum((gfibsval{side}(:,gpatsel)>obj.efieldthreshold),2);
        end
        % remove fibers that are not connected to enough VTAs/Efields or connected
        % to too many VTAs (connthreshold slider)
        gfibsval{side}(sumgfibsval<((obj.connthreshold/100)*length(gpatsel)),:)=0;
        % only in case of VTAs (given two-sample-t-test statistic) do we
        % need to also exclude if tract is connected to too many VTAs:
        if obj.statmetric==1
            gfibsval{side}(sumgfibsval>((1-(obj.connthreshold/100))*length(gpatsel)),:)=0;
        end
        
        switch obj.statmetric
            case 1 % t-tests
                % check if covariates exist:
                if exist('covars', 'var')
                    % they do:
                    vals{group,side}=nan(size(gfibsval{side},1),1);
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
                        vals{group,side}=ea_corrsignan(vals{group,side},ps',obj);
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
                vals{group,side}=nan(size(gfibsval{side},1),1);
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
                            outvals=partialcorr(invals,I(gpatsel,side),usecovars,'rows','pairwise','type',obj.corrtype); % generate optimality values on all but left out patients
                        else % no need to calc p-val here
                            [outvals,ps]=partialcorr(invals,I(gpatsel,side),usecovars,'rows','pairwise','type',obj.corrtype); % generate optimality values on all but left out patients
                            outvals=ea_corrsignan(outvals,ps,obj);
                        end
                    else
                        if conventionalcorr
                            if obj.showsignificantonly
                                [outvals,ps]=corr(invals,I(gpatsel,side),'rows','pairwise','type',obj.corrtype); % generate optimality values on all but left out patients
                                outvals=ea_corrsignan(outvals,ps,obj);
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
                                        [outvals,ps]=ea_bendcorr(invals,I(gpatsel,side)); % generate optimality values on all but left out patients
                                        outvals=ea_corrsignan(outvals,ps,obj);
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
                end
        end
        
        fibcell{group,side}=obj.results.(ea_conn2connid(obj.connectome)).fibcell{side}(~isnan(vals{group,side}));
        % Remove vals and fibers outside the thresholding range
        
        obj.stats.pos.available(side)=sum(vals{1,side}>0); % only collected for first group (positives)
        obj.stats.neg.available(side)=sum(vals{1,side}<0);
        usedidx{group,side}=find(~isnan(vals{group,side}));
        vals{group,side}=vals{group,side}(usedidx{group,side}); % final weights for surviving fibers
    end
    
    allvals = vertcat(vals{group,:});
    posvals = sort(allvals(allvals>0),'descend');
    negvals = sort(allvals(allvals<0),'ascend');
    
    for side=1:numel(gfibsval)
        if ~obj.posvisible || ~obj.showposamount(side) || isempty(posvals)
            posthresh = inf;
        else
            posrange = posvals(1) - posvals(end);
            posthresh = posvals(1) - obj.showposamount(side)/100 * posrange;
            
            if posrange == 0
                posthresh = posthresh - eps*10;
            end
        end
        
        if ~obj.negvisible || ~obj.shownegamount(side) || isempty(negvals)
            negthresh = -inf;
        else
            negrange = negvals(1) - negvals(end);
            negthresh = negvals(1) - obj.shownegamount(side)/100 * negrange;
            
            if negrange == 0
                negthresh = negthresh + eps*10;
            end
        end
        
        % Remove vals and fibers outside the thresholding range
        remove = logical(logical(vals{group,side}<posthresh) .* logical(vals{group,side}>negthresh));
        vals{group,side}(remove)=[];
        fibcell{group,side}(remove)=[];
        usedidx{group,side}(remove)=[];
    end
end



function vals=ea_corrsignan(vals,ps,obj)

nnanidx=~isnan(vals);
numtests=sum(nnanidx);

switch lower(obj.multcompstrategy)
    case 'fdr'
        pnnan=ps(nnanidx);
        [psort,idx]=sort(pnnan);
        pranks=zeros(length(psort),1);
        for rank=1:length(pranks)
           pranks(idx(rank))=rank; 
        end
        pnnan=pnnan.*numtests;
        pnnan=pnnan./pranks;
        ps(nnanidx)=pnnan;
    case 'bonferroni'
        ps(nnanidx)=ps(nnanidx).*numtests;
end
ps(~nnanidx)=1;
vals(ps>obj.alphalevel)=nan; % delete everything nonsignificant.














