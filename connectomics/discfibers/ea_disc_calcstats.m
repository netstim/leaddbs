function [vals]=ea_disc_calcstats(fibsval,I,obj,patsel)

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
        gpatsel=[gpatsel;gpatsel+length(obj.allpatients)];
    end
    
    for side=1:2
        
        % check connthreshold
        switch obj.statmetric
            case 1
                sumgfibsval=sum(gfibsval{side}(:,gpatsel),2);
            case 2
                sumgfibsval=sum((gfibsval{side}(:,gpatsel)>obj.efieldthreshold),2);
        end
        gfibsval{side}(sumgfibsval<((obj.connthreshold/100)*length(gpatsel)),:)=0;
        switch obj.statmetric
            case 1 % t-tests
                allvals=repmat(I(gpatsel,side),1,size(gfibsval{side}(:,gpatsel),1)); % improvement values (taken from Lead group file or specified in line 12).
                fibsimpval=allvals; % Make a copy to denote improvements of connected fibers
                fibsimpval(~logical(gfibsval{side}(:,gpatsel)))=nan; % Delete all unconnected values
                nfibsimpval=allvals; % Make a copy to denote improvements of unconnected fibers
                nfibsimpval(logical(gfibsval{side}(:,gpatsel)))=nan; % Delete all connected values
                [~,p,~,stats]=ttest2(fibsimpval,nfibsimpval); % Run two-sample t-test across connected / unconnected values
                vals{group,side}=stats.tstat;
                %vals{group,side}(p>0.5)=nan; % discard noisy fibers (optional or could be adapted)
            case 2 % spearmans correlations
                nangfibsval=gfibsval{side}(:,gpatsel);
                nangfibsval(nangfibsval==0)=nan; % only used in spearmans correlations
                vals{group,side}=corr(nangfibsval',I(gpatsel,side),'rows','pairwise','type','Spearman'); % generate optimality values on all but left out patients
        end
    end
end
