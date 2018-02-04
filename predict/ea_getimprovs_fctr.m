function I=ea_getimprovs_fctr(cdat,cmd,pts,fctr)

% returns improvement on factor supplied by cdat struct

sessions={'baseline','postop'};
for ses=1:length(sessions)
    factors.(sessions{ses})=nan(size(cdat.baseline,1),length(cdat.factornames));
    cdat.(sessions{ses})(isnan(cdat.(sessions{ses})))=0;
    
    for sc=1:length(cdat.scores)
        factors.(sessions{ses})(ismember(cdat.score,cdat.scores{sc}),:)=...
            cdat.(sessions{ses})(ismember(cdat.score,cdat.scores{sc}),:)*cdat.factors.(cdat.scores{sc});
    end
end

factors.baseline(factors.baseline==0)=nan;
switch cmd
    case 'raw'
        factors.improvement=(factors.baseline-factors.postop); %./factors.baseline;
        factors.improvement(isinf(factors.improvement))=nan;
        improvs=factors.improvement(:,fctr);        
        I=improvs(pts);
    case 'percent'
        factors.improvement=(factors.baseline-factors.postop)./factors.baseline;
        factors.improvement(isinf(factors.improvement))=nan;
        improvs=factors.improvement(:,fctr);
        I=improvs(pts);
    case 'baseline'
        improvs=factors.baseline(:,fctr);
        I=improvs(pts);
end