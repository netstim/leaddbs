function I=ea_getimprovs_fctr_smtp(cdat,cmd,pts,fctr,somatotopy)

% returns improvement on factor supplied by cdat struct
warning off
if ~exist('pts','var') || isempty(pts)
    pts=1:size(cdat.baseline,1); % all patients.
end
if ~exist('fctr','var')
    fctr=1;
end
if ~exist('somatotopy','var')
    somatotopy=1;
end
if ischar(fctr)
    [~,ix]=ismember(fctr,cdat.factornames);
    if ix % found match
        fctr=ix;
    else
        fctrs=strsplit(fctr,'*'); % create new combined factor on the fly (e.g. could be 'anxiosomatic*mood')
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


if ischar(somatotopy)
    [~,ix]=ismember(somatotopy,cdat.somatotopynames);
    if ix % found match
        somatotopy=ix;
    else
        somatotopies=strsplit(somatotopy,'*'); % create new combined somatotopy on the fly (e.g. could be 'leftbody*upper_extremity')
        cdat.somatotopynames{end+1}='custom';
        for sc=1:length(cdat.scores)
            cdat.somatotopies.(cdat.scores{sc})(:,end+1)=ones(size(cdat.somatotopies.(cdat.scores{sc}),1),1);
            for fc=1:length(somatotopies)
                [~,ix]=ismember(somatotopies{fc},cdat.factornames);
                if ix
                cdat.somatotopies.(cdat.scores{sc})(:,end)=cdat.somatotopies.(cdat.scores{sc})(:,end).*...
                    cdat.somatotopies.(cdat.scores{sc})(:,ix);
                end
            end
        end     
        somatotopy=length(cdat.somatotopynames);
    end
end


% factor and somatotopy indices selected. now go select that subscore:

sessions={'baseline','postop'};
for ses=1:length(sessions)
    factors.(sessions{ses})=nan(size(cdat.baseline,1),1);
    cdat.(sessions{ses})(isnan(cdat.(sessions{ses})))=0;
    
    for sc=1:length(cdat.scores)
        factors.(sessions{ses})(ismember(cdat.score,cdat.scores{sc}),:)=...
            cdat.(sessions{ses})(ismember(cdat.score,cdat.scores{sc}),:)*...
            (cdat.factors.(cdat.scores{sc})(:,fctr).*...
            cdat.somatotopies.(cdat.scores{sc})(:,somatotopy));
    end
end

factors.baseline(factors.baseline<(max(factors.baseline)*0.1))=nan; % assure at least 10% symptoms are present
switch cmd
    case {'raw','absolute','abs'}
        factors.improvement=(factors.baseline-factors.postop); %./factors.baseline;
        factors.improvement(isinf(factors.improvement))=nan;
        improvs=factors.improvement;        
        I=improvs(pts);
    case {'percent','perc','rel','relative'}
        factors.improvement=(factors.baseline-factors.postop)./factors.baseline;
        factors.improvement(isinf(factors.improvement))=nan;
        improvs=factors.improvement;
        I=improvs(pts);
    case 'baseline'
        improvs=factors.baseline;
        I=improvs(pts);
    case 'cleaned' % cleaned from baseline by regressing out baseline.
        factors.improvement=(factors.baseline-factors.postop); %./factors.baseline;
        factors.improvement(isinf(factors.improvement))=nan;
        improvs=ea_resid(factors.baseline,factors.improvement);     
        I=improvs(pts);
        
end
warning on