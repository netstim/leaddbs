function result=ea_mahalanobis_cov(x,val)

% This function is the same as ea_mahalanobis but only optimized for a
% specific input. It is called from ea_mcdcov
%
% Is is equal to ea_mahalanobis(x,[0; 0],'cov',val) with size(val) = [2 1]

result=sum(x*pinv(diag(val)).*x,2)'; 
