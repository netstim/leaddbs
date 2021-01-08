function [r]=ea_skipped_correlation(x,y,type)

% performs a robust correlation using pearson/spearman correlation on
% data cleaned up for bivariate outliers - that is after finding the
% central point in the distribution using the mid covariance determinant,
% orthogonal distances are computed to this point, and any data outside the
% bound defined by the idealf estimator of the interquartile range is removed.
% 
% FORMAT:
%          [r,t,h] = skipped_correlation(X,Y);
%          [r,t,h] = skipped_correlation(X,Y,fig_flag);
%          [r,t,h,outid,hboot,CI] = skipped_correlation(X,Y,fig_flag);
%
% INPUTS:  X and Y are 2 vectors or matrices, in the latter case,
%          correlations are computed column-wise 
%          fig_flag (1/0) indicates to plot the data or not
%
% OUTPUTS:
%          r is the pearson/spearman correlation 
%          t is the T value associated to the skipped correlation
%          h is the hypothesis of no association at alpha = 5% 
%          outid is the index of bivariate outliers
% 
%          optional:
%
%          hboot 1/0 declares the test significant based on CI (h depends on t)
%          CI is the robust confidence interval computed by bootstrapping the 
%          cleaned-up data set and taking the .95 centile values
%
% This code rely on the mid covariance determinant as implemented in LIBRA
% - Verboven, S., Hubert, M. (2005), LIBRA: a MATLAB Library for Robust Analysis,
% Chemometrics and Intelligent Laboratory Systems, 75, 127-136.
% - Rousseeuw, P.J. (1984), "Least Median of Squares Regression,"
% Journal of the American Statistical Association, Vol. 79, pp. 871-881.
% 
% The quantile of observations whose covariance is minimized is 
% floor((n+size(X,2)*2+1)/2)),
% i.e. ((number of observations + number of variables*2)+1) / 2, 
% thus for a correlation this is floor(n/2 + 5/2).
%
% See also MCDCOV, IDEALF.

% Cyril Pernet & Guillaume Rousselet, v1 - April 2012
% ---------------------------------------------------
%  Copyright (C) Corr_toolbox 2012

% transpose if x or y are not in column
if size(x,1) == 1 && size(x,2) > 1; x = x'; end
if size(y,1) == 1 && size(y,2) > 1; y = y'; end

if ~exist('type','var')
    type='Pearson';
end

% now if x is a vector and we test multiple y (or the other way around) one
% has to adjust for this
if size(x,2) == 1 && size(y,2) > 1
    x = repmat(x,1,size(y,2));
elseif size(y,2) == 1 && size(x,2) > 1
    y = repmat(y,1,size(x,2));
end

[n,p] = size(x);
if size(x) ~= size(y)
    error('x and y are of different sizes')
elseif n < 10
    error('robust effects can''t be computed with less than 10 observations')
end

gval = sqrt(chi2inv(0.975,2)); % in fact depends on size(X,2) but here always = 2

%% compute
for column = 1:p
    
    % get the centre of the bivariate distributions
    X = [x(:,column) y(:,column)];
    X = X(~logical(sum(isnan(X),2)),:); % Remove NAN
    N = size(X,1); % New row number after removing NAN
    if N < 10 % Check again after removing NAN
        error('robust effects can''t be computed with less than 10 observations')
    end

    result=ea_mcdcov(X,'cor',1,'plots',0,'h',floor((N+size(X,2)*2+1)/2));
    center = result.center;

    % orthogonal projection to the lines joining the center
    % followed by outlier detection using mad median rule
    
    vec=1:N;
    for i=1:N % for each row
        dis=NaN(N,1);
        B = (X(i,:)-center)';
        BB = B.^2;
        bot = sum(BB);
        if bot~=0
            for j=1:N
                A = (X(j,:)-center)';
                dis(j)= norm(A'*B/bot.*B); 
            end
            % MAD median rule
            %[outliers,value] = madmedianrule(dis,2);
            %record{i} = dis > (median(dis)+gval.*value);
            % IQR rule
            [ql,qu]=ea_idealf(dis);
            record{i} = (dis > median(dis)+gval.*(qu-ql)) ; % + (dis < median(dis)-gval.*(qu-ql));
        end
    end
    
    try
        flag = nan(N,1);
        flag = sum(cell2mat(record),2); % if any point is flagged
        
    catch ME  % this can happen to have an empty cell so loop
        flag = nan(N,size(record,2));
        index = 1;
        for s=1:size(record,2)
          if ~isempty(record{s})
              flag(:,index) = record{s};
              index = index+1;
          end
        end
        flag(:,index:end) = [];
        flag = sum(flag,2);
    end
            
    if sum(flag)==0
        outid{column}=[];
    else
        flag=(flag>=1);
        outid{column}=vec(flag);
    end
    keep=vec(~flag);
    
    %% Pearson/Spearman correlation
    a = X(keep,1);
    b = X(keep,2);

    switch lower(type)
        case 'pearson'
            r(column) = sum(detrend(a,'constant').*detrend(b,'constant')) ./ ...
                (sum(detrend(a,'constant').^2).*sum(detrend(b,'constant').^2)).^(1/2);
        case 'spearman'
            xrank = tiedrank(a,0); yrank = tiedrank(b,0);
            r(column) = sum(detrend(xrank,'constant').*detrend(yrank,'constant')) ./ ...
                (sum(detrend(xrank,'constant').^2).*sum(detrend(yrank,'constant').^2)).^(1/2);
    end

    %disp(num2str(r(column)));
end
