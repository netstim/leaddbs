% fweAdjust.m
%
%        $Id$
%      usage: adjustedP = fweAdjust(p,params,)
%         by: julien besle, 
%       date: 17/01/2011
%    purpose: adjusts p-values for various familywise error control procedure
%            (Bonferroni, Holm, Hochberg and Hommel)
%            for adjustment versions of each test, see Wright (1992) BIOMETRICS 48, 1005-1013 
%            In addition, the adjustment can be made more powerful by providing an estimate
%            of the proportion of true null hypotheses and its corresponding unadjusted threshold 
%            (see estimateNumberTrueH0.m  and see Sarkar, 2009, submitted)
%

function adjustedPdata = fweAdjust(p, params, numberTrueH0, lambda)


isNotNan = ~isnan(p);
sizePdata = size(p);
p = p(isNotNan);
[p,sortingIndex] = sort(p);
numberH0 = length(p);

if ieNotDefined('numberTrueH0') || ieNotDefined('lambda')
  numberTrueH0 = numberH0;
  lambda = 0;
end


switch(params.fweMethod)
  case {'Single-step (Bonferroni)','Adaptive Single-step'} %Bonferroni (adjustment version)
    adjustedP = p*numberTrueH0; %no need to use sorted data, but makes the code simpler
    
  case {'Step-down (Holm)','Adaptive Step-down'} %Holm's procedure (adjustment version)
    %Holm (1979) Scandinavian Journal of Statistics, Vol. 6, No. 2 , pp. 65-70
    adjustedP = p.*(1+lambda).*min((numberH0:-1:1),numberTrueH0)';
    %enforce monotonicity by taking the max value between a p-value and the previous (smaller) one
    for i=2:numberH0
     adjustedP(i) = max(adjustedP(i),adjustedP(i-1));
    end
    
  case 'Step-up (Hochberg)' %Hochberg's procedure (adjustment version)
    %Hochberg (1988) Biometrika, 75, 4, pp. 800-2
    adjustedP = p.*(1+lambda).*min((numberH0:-1:1),numberTrueH0)';
    %enforce monotonicity by taking the min between a p-value and the next (larger) one
    for i=numberH0-1:-1:1
     adjustedP(i) = min(adjustedP(i),adjustedP(i+1));
    end
    
  case 'Hommel' 
    % Hommel (1988) Biometrika (1988), 75, 2, pp. 383-6
%%%    % p-adjustment algorithm provided by Wright (1992)
% % % %     adjustedP=p;
% % % %     for m=numberH0:-1:2
% % % %       cMin = min(m*p(numberH0-m+1:numberH0)./(m+(numberH0-m+1:numberH0)-numberH0));
% % % %       adjustedP(numberH0-m+1:numberH0) = max(adjustedP(numberH0-m+1:numberH0),cMin);
% % % %       adjustedP(1:numberH0-m) = max(adjustedP(1:numberH0-m),min(cMin,m*p(1:numberH0-m)));
% % % %     end
    
    %I think this version is clearer:
    adjustedP=p;
    for i=1:numberH0-1
      cMin = min((numberH0+1-i)*p(i:numberH0)./(1:numberH0+1-i)');
      adjustedP(i:numberH0) = max(adjustedP(i:numberH0),cMin);
      adjustedP(1:i-1) = max(adjustedP(1:i-1),min(cMin,(numberH0+1-i)*p(1:i-1)));
    end

    
end

adjustedP = min(adjustedP,1); %bound with 1

adjustedP(sortingIndex) = adjustedP;
adjustedPdata = NaN(sizePdata);
adjustedPdata(isNotNan) = adjustedP;
