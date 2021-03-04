function [r,p] = ea_bendcorr(X,Y)

% Computes the percentage bend correlation along with the bootstrap CI
%
% FORMAT:  [r,t,p] = bendcorr(X,Y)
%          [r,t,p,hboot,CI,H,pH] = bendcorr(X,Y,fig_flag,beta)
%
% INPUTS:  X and Y are 2 vectors or matrices. In the latter case,
%          correlations are computed column-wise. 
%          fig_flag indicates to plot (1 - default) the data or not (0)
%          beta represents the amount of trimming: 0 <= beta <= 0.5
%          (beta is also called the bending constant for omega - default = 0.2)
%
% OUTPUTS: r is the percentage bend correlation
%          t is the associated t value
%          pval is the corresponding p value
%          hboot 1/0 declares the test significant based on CI
%          CI is the percentile bootstrap confidence interval
%          H is the measure of association between all pairs
%          pH is the p value for an omnibus test of independence between all pairs 
%
% The percentage bend correlation is a robust method that protects against
% outliers among the marginal distributions.

% Cyril Pernet and Guillaume Rousselet 26-01-2011
% Reformatted for Corr_toolbox 02--7-2012 
% ----------------------------------------------
%  Copyright (C) Corr_toolbox 2012

%% data check

if nargin<2
    error('two input vectors requested')
elseif nargin>4
    eror('too many inputs')
end

% if X a vector and Y a matrix, 
% repmat X to perform multiple tests on Y (or the other around)
if size(X,1) == 1 && size(X,2) > 1; X = X'; end
if size(Y,1) == 1 && size(Y,2) > 1; Y = Y'; end

if size(X,2) == 1 && size(Y,2) > 1
    X = repmat(X,1,size(Y,2));
elseif size(Y,2) == 1 && size(X,2) > 1
    Y = repmat(Y,1,size(X,2));
end

if sum(size(X)~=size(Y)) ~= 0
    error('X and Y must have the same size')
end

%% parameters
level = 5/100;
beta=0.5;



% remove NaNs
% -----------
X = [X Y];
X(find(sum(isnan(X),2)),:) = [];
n = length(X);

%% compute
% --------

    [r,p] = bend_compute(X,beta);


r=r';
p=p';
 


end

function [r,p] = bend_compute(X,beta)

H= []; pH = [];

%% Medians and absolute deviation from the medians
% ---------------------------------------------
M = repmat(median(X),size(X,1),1);
W = sort(abs(X-M),1); 
 
% limits
% -------
m = floor((1-beta)*size(X,1));
omega = W(m,:); 

%% Compute the correlation
% ------------------------
P = (X-M)./ repmat(omega,size(X,1),1); 
P(isnan(P)) = 0; P(isinf(P)) = 0; % correct if omega = 0
comb = [(1:size(X,2)/2)',((1:size(X,2)/2)+size(X,2)/2)']; % all pairs of columns
r = NaN(size(comb,1),1); 
t = r; p = t;
for j = 1:size(comb,1)
    
    % column 1
    psi = P(:,comb(j,1)); 
    i1 = length(psi(psi<-1)); 
    i2 = length(psi(psi>1)); 
    sx = X(:,comb(j,1)); 
    sx(psi<(-1)) = 0; 
    sx(psi>1) = 0; 
    pbos = (sum(sx)+ omega(comb(j,1))*(i2-i1)) / (size(X,1)-i1-i2); 
    a = (X(:,comb(j,1))-pbos)./repmat(omega(comb(j,1)),size(X,1),1); 
        
    % column 2
    psi = P(:,comb(j,2));
    i1 = length(psi(psi<-1));
    i2 = length(psi(psi>1));
    sx = X(:,comb(j,2));
    sx(psi<(-1)) = 0;
    sx(psi>1) = 0;
    pbos = (sum(sx)+ omega(comb(j,2))*(i2-i1)) / (size(X,1)-i1-i2);
    b = (X(:,comb(j,2))-pbos)./repmat(omega(comb(j,2)),size(X,1),1);
    
     % bend
     a(a<=-1) = -1; a(a>=1) = 1; 
     b(b<=-1) = -1; b(b>=1) = 1; 
     
     % get r & p
     r(j) = sum(a.*b)/sqrt(sum(a.^2)*sum(b.^2));
     t(j) = r(j)*sqrt((size(X,1) - 2)/(1 - r(j).^2));
     p(j) = 2*(1 - tcdf(abs(t(j)),size(X,1)-2));

end



end