function [vals,fibcell,usedidx] = ea_discfibers_calcstats(obj,patsel,Iperm)

if ~exist('Iperm','var')
    I=obj.responsevar;
else % used in permutation based statistics - in this case the real improvement can be substituted with permuted variables.
    I=Iperm;
end

if obj.multresponsevarneg % flag to multiply response var by -1
    I=-I;
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

    for side=1:2
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
                    [~,~,~,stats]=ttest2(fibsimpval',nfibsimpval'); % Run two-sample t-test across connected / unconnected values
                    vals{group,side}=stats.tstat';
                    %vals{group,side}(p>0.5)=nan; % discard noisy fibers (optional or could be adapted)
                end
            case 2 % spearmans correlations
                 if exist('covars', 'var')
                    usecovars=[];
                    for cv=1:length(covars)
                        thiscv=covars{cv}(gpatsel,:);
                        if (size(thiscv,2)==2)
                            thiscv=thiscv(:,side);
                        end
                        usecovars=[usecovars,thiscv];
                    end
                    vals{group,side}=partialcorr(gfibsval{side}(:,gpatsel)',I(gpatsel,side),usecovars,'rows','pairwise','type','Spearman'); % generate optimality values on all but left out patients
                else
                    vals{group,side}=corr(gfibsval{side}(:,gpatsel)',I(gpatsel,side),'rows','pairwise','type','Spearman'); % generate optimality values on all but left out patients
                end
        end

        fibcell{group,side}=obj.results.(ea_conn2connid(obj.connectome)).fibcell{side}(~isnan(vals{group,side}));
        % Remove vals and fibers outside the thresholding range

        obj.stats.pos.available(side)=sum(vals{1,side}>0); % only collected for first group (positives)
        obj.stats.neg.available(side)=sum(vals{1,side}<0);
        usedidx{group,side}=find(~isnan(vals{group,side}));
        vals{group,side}=vals{group,side}(usedidx{group,side}); % final weights for surviving fibers
    end

    allvals = [vals{group,1};vals{group,2}];
    posvals = sort(allvals(allvals>0),'descend');
    posrange = posvals(1) - posvals(end);
    negvals = sort(allvals(allvals<0),'ascend');
    negrange = negvals(1) - negvals(end);

    for side=1:2
        if ~obj.posvisible
            posthresh = posvals(1);
        else
            try
                posthresh = posvals(1) - obj.showposamount(side)/100 * posrange;
            catch
                posthresh = posvals(1);
            end
        end
        posthresh = posthresh + eps;

        if ~obj.negvisible
            negthresh = negvals(1);
        else
            try
                negthresh = negvals(1) - obj.shownegamount(side)/100 * negrange;
            catch
                negthresh = negvals(1);
            end
        end
        negthresh = negthresh - eps;

        % Remove vals and fibers outside the thresholding range
        remove = logical(logical(vals{group,side}<posthresh) .* logical(vals{group,side}>negthresh));
        vals{group,side}(remove)=[];
        fibcell{group,side}(remove)=[];
        usedidx{group,side}(remove)=[];
    end
end
