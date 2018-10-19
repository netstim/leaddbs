function res=ea_resid(X,y)
% gives out residuals directly.


%X(isnan(X))=0;
%y(isnan(y))=0;
[b,dev,stats]=glmfit(X,y);
res=stats.resid;