function ypred=ea_lo_prediction(X,y,modix,predix,analys,corrtype)

if ~exist('corrtype','var')
corrtype='Spearman';
end
if ~exist('analys','var')
    analys={'R','W','C'}; % do all three types
end

N=length(y);


% build model on modix
if ismember('R',analys) || ismember('C',analys)
    R=corr(X(modix,:),y(modix),'rows','pairwise','type',corrtype)'; % R-Map
end
if ismember('W',analys) || ismember('C',analys)
    sX=nansum(X(modix,:).*repmat(y(modix),1,size(X,2),1)); % Weighted average Map
end
if ismember('C',analys)
    CsX=sX;
    for ix=1:size(X,2)
        if (CsX(ix)>0 && R(ix)<0) || (CsX(ix)<0 && R(ix)>0)
            CsX(ix)=nan;
        end
    end
end


for leo=predix    
    iterpt=1:N;
    iterpt(leo)=[];
    for pt=1:N % iterate through all here.
        Rr(pt)=corr(X(pt,:),R(:),'rows','pairwise','type',corrtype);
        WavgRr(pt)=corr(nii.img(:),sX(:),'rows','pairwise','type',corrtype);
        CsXr(pt)=corr(nii.img(:),CsX(:),'rows','pairwise','type',corrtype);
    end
    
    % fit model to predict improvs based on other patients:
    if ismember('R',analys)
        b=glmfit(Rr(iterpt),y(iterpt));
        ypred.R(leo)=ea_addone(Rr(leo))*b;
    end
    if ismember('W',analys)
        b=glmfit(WavgRr(iterpt),y(iterpt));
        ypred.W(leo)=ea_addone(WavgRr(leo))*b;
    end
    if ismember('C',analys)
        b=glmfit(CsXr(iterpt),y(iterpt));
        ypred.C(leo)=ea_addone(CsXr(leo))*b;
    end
end





    