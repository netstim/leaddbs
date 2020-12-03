function [vals,usedidx] = ea_networkmapping_calcstats(obj,patsel,Iperm)

if ~exist('Iperm','var')
    I=obj.responsevar;
else % used in permutation based statistics - in this case the real improvement can be substituted with permuted variables.
    I=Iperm;
end

AllX = (obj.results.(ea_conn2connid(obj.connectome)).connval);

% quickly recalc stats:
if ~exist('patsel','var') % patsel can be supplied directly (in this case, obj.patientselection is ignored), e.g. for cross-validations.
    patsel=obj.patientselection;
end

if size(I,2)==2 % 1 entry per patient, not per electrode
    I=mean(I,2); % both sides the same;
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
    if dogroups
        groupspt=find(obj.M.patient.group==group);
        gpatsel=groupspt(ismember(groupspt,patsel));
    else
        gpatsel=patsel;
    end
    if obj.mirrorsides
        gpatsel=[gpatsel,gpatsel+length(obj.allpatients)];
    end
    

        switch obj.statmetric
            case 'Correlations (R-map)' 
                % check if covariates exist:
                if exist('covars', 'var')
                    % they do:
                    keyboard % need to add in
                else
                    % no covariates exist:
                    if obj.showsignificantonly
                        [vals{group},ps]=corr(I(gpatsel),(AllX(gpatsel,:)),'rows','pairwise','type',obj.corrtype); % improvement values (taken from Lead group file or specified in line 12).
                        vals{group}=ea_corrsignan(vals{group},ps',obj);
                    else
                        [vals{group}]=corr(I(gpatsel),(AllX(gpatsel,:)),'rows','pairwise','type',obj.corrtype); % improvement values (taken from Lead group file or specified in line 12).
                    end
                end
            case 'Weighted ave..' % check
                keyboard
        end
                
        obj.stats.pos.available=sum(vals{1}>0); % only collected for first group (positives)
        obj.stats.neg.available=sum(vals{1}<0);
        usedidx{group}=find(~isnan(vals{group}));
        vals{group}=vals{group}(usedidx{group}); % final weights for surviving fibers    
%         allvals = vertcat(vals{group,:});
%         posvals = sort(allvals(allvals>0),'descend');
%         negvals = sort(allvals(allvals<0),'ascend');
%         
%         if ~obj.posvisible || ~obj.showposamount || isempty(posvals)
%             posthresh = inf;
%         else
%             posrange = posvals(1) - posvals(end);
%             posthresh = posvals(1) - obj.showposamount/100 * posrange;
%             
%             if posrange == 0
%                 posthresh = posthresh - eps*10;
%             end
%         end
%         
%         if ~obj.negvisible || ~obj.shownegamount || isempty(negvals)
%             negthresh = -inf;
%         else
%             negrange = negvals(1) - negvals(end);
%             negthresh = negvals(1) - obj.shownegamount/100 * negrange;
%             
%             if negrange == 0
%                 negthresh = negthresh + eps*10;
%             end
%         end
%         
%         % Remove vals and fibers outside the thresholding range
%         remove = logical(logical(vals{group}<posthresh) .* logical(vals{group}>negthresh));
%         vals{group}(remove)=[];
%         usedidx{group}(remove)=[];
    
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














