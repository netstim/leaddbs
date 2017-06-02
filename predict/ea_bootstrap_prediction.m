function [ypredR,ypredW,ypredC]=ea_bootstrap_prediction(X,y,modix,predix,corrtype,regtype)

if ~exist('corrtype','var')
    corrtype='Spearman';
end
if ~exist('regtype','var')
    regtype='nonparametric';
end

N=length(y);


% build model on modix
R=corr(X(modix,:),y(modix),'rows','pairwise','type',corrtype)'; % R-Map
if nargout>1
    sX=nansum(X(modix,:).*repmat(y(modix),1,size(X,2),1)); % Weighted average Map
end
if nargout>2
    CsX=sX;
    CsX((CsX.*R)<0)=nan;
end

% fit model to predict improvs based on model patients:
Rr=corr(X',R','rows','pairwise','type',corrtype);
switch lower(regtype)
    case 'parametric'
        bR=glmfit(Rr(modix),y(modix));
    case 'nonparametric'
        ksrR=ea_ksr(Rr(modix),y(modix));
end
if nargout>1
    WavgRr=corr(X',sX','rows','pairwise','type',corrtype);
    switch lower(regtype)
        case 'parametric'
            bW=glmfit(WavgRr(modix),y(modix));
        case 'nonparametric'
            ksrW=ea_ksr(WavgRr(modix),y(modix));
    end
end
if nargout>2
    CsXr=corr(X',CsX','rows','pairwise','type',corrtype);
    switch lower(regtype)
        case 'parametric'
            bC=glmfit(CsXr(modix),y(modix));
        case 'nonparametric'
            ksrC=ea_ksr(CsXr(modix),y(modix));
    end
end


%cnt=1;
switch lower(regtype)
    case 'parametric'
        ypredR=ea_addone(Rr(predix))*bR;
        
        if nargout>1
            ypredW=ea_addone(WavgRr(predix))*bW;
        end
        if nargout>2
            ypredC=ea_addone(CsXr(predix))*bC;
        end
    case 'nonparametric'
        ypredR=interp1(ksrR.x,ksrR.f,Rr(predix));
        
        if nargout>1
            ypredW=interp1(ksrW.x,ksrW.f,WavgRr(predix));
        end
        if nargout>2
            ypredC=interp1(ksrC.x,ksrC.f,CsXr(predix));
        end
end


%RRR=corr([y(predix),ypredR,ypredW,ypredC]);
%disp(['Bootstrap prediction yields R = ',num2str(RRR(1,2)),' (R-maps), ',num2str(RRR(1,3)),' (Weighted average maps), ',num2str(RRR(1,4)),' (Combined maps).']);





