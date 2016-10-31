% estimateNumberTrueH0.m
%
%        $Id$
%      usage: adjustedP = estimateNumberTrueH0(p,params,)
%         by: julien besle, 
%       date: 17/01/2011
%    purpose: estimates the proportion of null hypotheses (see Sarkar, 2009, submitted)
%

function [numberTrueH0,lambda] = estimateNumberTrueH0(p, params)

isNotNan = ~isnan(p);
%sizePdata = size(p);
p = p(isNotNan);
[p,sortingIndex] = sort(p);
numberH0 = length(p);


switch(params.trueNullsEstimationMethod)
  case 'Raw P-values Cut-off'
    %Storey et al. (2004) Sarkar, 2009, submitted
    lambda = params.trueNullsEstimationThreshold;
    numberTrueH0 = (numberH0 - nnz(p<lambda)+1)/(1-lambda);
    numberTrueH0 = min(numberTrueH0,numberH0);

  case 'Lowest Slope'
    %Hochberg and Benjamini (1990), Benjamini and Hochberg, (2000)
    %compute the estimated number of true H0 when considering the number of rejections at all possible levels
    numberTrueH0 = (numberH0+1-(1:numberH0)')./(1-p);
    %find the first estimated number that is greater than the previous one
    k = find(diff(numberTrueH0)>0,1,'first');
    numberTrueH0 = min(numberTrueH0(k+1),numberH0);
    lambda = p(k); %this is the raw p-value cut-off that would correspond 
    %to the estimated number of true H0 using Storey's method

  case 'Least Squares'
    %Hsueh et al. (2003) Journal of biopharmaceutical statistics 13(4), p675–689
    %compute the power at each possible raw p-value cutoff (effectively the actual p-values provided)
    beta = cumsum(1:numberH0)'./(1:numberH0)'.^2;
    %the expected number of rejection for each possible cutoff/raw p-value is
    % numberH0*p(i) + (numberH0 - numberTrueH0)*beta(i) for i=1:numberH0
    % find the value of numberTrueH0 that minimizes the sum of square difference 
    % between this expected number and the actual number of rejections for each cutoff (1:numberH0)
    %that is sum( (i - numberH0*p(i) - (numberH0 - numberTrueH0)*beta(i))^2 )
    %this is given by
    x = p-beta;
    y = (1:numberH0)' - numberH0*beta;
    numberTrueH0 = sum(x.*y)/sum(x.^2);
    numberTrueH0 = min(numberTrueH0,numberH0);
    
    %for lambda, find the raw p-value cut-off that would correspond 
    %to the estimated number of true H0 using Storey's method
    %the problem is there are several possible values
    %so let's just take the most conservative
    lambda = p(find((numberH0+1-(1:numberH0)')./(1-p)<numberTrueH0,1,'first'));
    
  otherwise
    numberTrueH0 = numberH0;
    lambda =0;
end
