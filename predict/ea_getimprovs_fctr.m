function I=ea_getimprovs_fctr(cdat,cmd,pts,fctr)

% returns improvement on factor supplied by cdat struct
warning off
if ~exist('pts','var') || isempty(pts)
    pts=1:size(cdat.baseline,1); % all patients.
end
if ~exist('fctr','var')
    fctr=1;
end
if ischar(fctr)
    [~,ix]=ismember(fctr,cdat.factornames);
    if ix % found match
        fctr=ix;
    else
        fctrs=strsplit(fctr,'*');
        cdat.factornames{end+1}='custom';
        for sc=1:length(cdat.scores)
            cdat.factors.(cdat.scores{sc})(:,end+1)=ones(size(cdat.factors.(cdat.scores{sc}),1),1);
            for fc=1:length(fctrs)
                [~,ix]=ismember(fctrs{fc},cdat.factornames);
                if ix
                cdat.factors.(cdat.scores{sc})(:,end)=cdat.factors.(cdat.scores{sc})(:,end).*...
                    cdat.factors.(cdat.scores{sc})(:,ix);
                end
            end
        end     
        fctr=length(cdat.factornames);
    end
end
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
    case {'raw','absolute','abs'}
        factors.improvement=(factors.baseline-factors.postop); %./factors.baseline;
        factors.improvement(isinf(factors.improvement))=nan;
        improvs=factors.improvement(:,fctr);        
        I=improvs(pts);
    case {'percent','perc','rel','relative'}
        factors.improvement=(factors.baseline-factors.postop)./factors.baseline;
        factors.improvement(isinf(factors.improvement))=nan;
        improvs=factors.improvement(:,fctr);
        I=improvs(pts);
    case 'baseline'
        improvs=factors.baseline(:,fctr);
        I=improvs(pts);
    case 'cleaned' % cleaned from baseline by regressing out baseline.
        factors.improvement=(factors.baseline-factors.postop); %./factors.baseline;
        factors.improvement(isinf(factors.improvement))=nan;
        improvs=ea_resid(factors.baseline,factors.improvement(:,fctr));     
        I=improvs(pts);
        
end
warning on