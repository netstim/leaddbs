function [vals]=ea_disc_calcstats(obj,patsel,Iperm)

if ~exist('Iperm','var')
I=obj.responsevar;
else % used in permutation based statistics - in this case the real improvement can be substituted with permuted variables.
    I=Iperm;
end
fibsval=obj.results.(ea_conn2connid(obj.connectome)).(ea_method2methodid(obj)).fibsval;

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
                    ea_dispercent(0,['Side ',num2str(side),': Calculating T-values (taking covariates into account)']);
                    nixfib=find(any(gfibsval{side}(:,gpatsel)').*~all(gfibsval{side}(:,gpatsel)'));
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
    end
end
