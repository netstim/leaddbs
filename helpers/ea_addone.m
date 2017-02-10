function X=ea_addone(X)
% simple function to add intercept term (e.g. for regression)
X=[ones(size(X,1),1),X];