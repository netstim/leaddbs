function [convertedStatistic, fdrAdjustedStatistic, fweAdjustedStatistic] = transformStatistic(p, outputPrecision, params)
%[convertedStatistic, fdrAdjustedStatistic, fweAdjustedStatistic] = transformStatistic(p, outputPrecision, params)

convertedStatistic = convertStatistic(p, params.testOutput, outputPrecision);

if params.fdrAdjustment
  fdrAdjustedP = p;
  for iTest = 1:size(p,4)
    fdrAdjustedP(:,:,:,iTest) = fdrAdjust(p(:,:,:,iTest),params);
  end
  fdrAdjustedStatistic = convertStatistic(fdrAdjustedP, params.testOutput, outputPrecision);
else
  fdrAdjustedStatistic = [];
end

if params.fweAdjustment
  fweAdjustedP = p;
  for iTest = 1:size(p,4)
    if ismember(params.fweMethod,{'Adaptive Step-down','Adaptive Single-step'})
      [numberTrueH0,lambda] = estimateNumberTrueH0(p(:,:,:,iTest),params);
    else
      numberTrueH0 = nnz(~isnan(p(:,:,:,iTest)));
      lambda = 0;
    end
    fweAdjustedP(:,:,:,iTest) = fweAdjust(p(:,:,:,iTest),params,numberTrueH0,lambda);
  end
  fweAdjustedStatistic = convertStatistic(fweAdjustedP, params.testOutput, outputPrecision);
else
  fweAdjustedStatistic = [];
end

function p = convertStatistic(p, outputStatistic, outputPrecision)

switch(outputStatistic)
  case 'Z value'
    %convert P value to Z value, 
    p = norminv(1-p);  
    p(p<0) = 0; % we're not interested in negative Z value
    %if there was no round-off error from cdf, we could do the following:
    %Z = max(-norminv(p),0);  %because the left side of norminv seems to be less sensitive to round-off errors,
    %we get -norminv(x) instead of norminv(1-x). also we'renot interested in negative Z value
  case '-log10(P) value'
    %convert P to -log10(P) 
    p = -log10(p);
    
end
p = eval([outputPrecision '(p)']);

function adjustedPdata = fdrAdjust(pData,params)

isNotNan = ~isnan(pData);
sizePdata = size(pData);
pData = pData(isNotNan);
[pData,sortingIndex] = sort(pData);
numberH0 = length(pData); 

switch(params.fdrAssumption)
  case 'Independence/Positive dependence'
    dependenceCorrectionConstant =1;
  case 'None'
    dependenceCorrectionConstant = sum(1./(1:numberH0));
    %approximate value: 
    %dependenceCorrectionConstant = log(numberH0)+0.57721566490153 %Euler constant
end

%estimate number of True Null hypotheses
% ref: Benjamini et al. 2006, Biometrika, 93(3) p491
switch(params.fdrMethod)
  case 'Step-up'
      adjustedP = linearStepUpFdr(pData,dependenceCorrectionConstant,1); 
      
  case 'Adaptive Step-up' %definition 3 in Benjamini et al. 2006
    %estimate the number of true H0 using one of the methods specified by params.trueNullsEstimationMethod
    numberTrueH0 = min(ceil(estimateNumberTrueH0(pData,params)),numberH0);
    %recompute the adjusted p-value using the estimated ratio of true H0
    adjustedP = linearStepUpFdr(pData,dependenceCorrectionConstant,numberTrueH0/numberH0);
    
  case 'Two-stage Step-up'
    %we will estimate the number of true H0 using definition 6 in Benjamini et al. 2006
    % first adjust the p values with q'=(q+1)
    adjustedP = linearStepUpFdr(pData,dependenceCorrectionConstant,1);
    numberFalseH0 = length(find(adjustedP<(params.trueNullsEstimationThreshold/(params.trueNullsEstimationThreshold+1))));
    if numberFalseH0 %if any voxel is significant at the chosen level usually .05/(.05+1)
      numberH0 = length(adjustedP); 
      %recompute the adjusted p-value using the estimated number of true H0
      adjustedP = linearStepUpFdr(pData,dependenceCorrectionConstant,(numberH0 - numberFalseH0)/numberH0);
    end
    
  case 'Multiple-stage Step-up'
    %we will estimate the number of true H0 using definition 7 in Benjamini et al. 2006
    %but they only give the definition for p-value correction with a fixed q
    q = params.trueNullsEstimationThreshold;
    %step-up version with L>=J (I think it's a step-up procedure within a step-down procedure)
    %k=max{i : for all j<=i there exists l>=j so that p(l)<=ql/{m+1?j(1?q)}}.
    i=0;          %
    l = numberH0; %these are just arbitrary values to pass the first while test
    while i<numberH0 && ~isempty(l)
      i = i+1;
      l = find(pData(i:numberH0)<=q*(i:numberH0)'./(numberH0+1-i*(1-q)),1,'last');
    end
    k=i;
    correctedP = pData;
    correctedP(k:end) = 1;
    % adjustement procedure seems a bit complicated... so I'll leave it for now
    % I don't think the results would be much different from the next, simpler, step-down procedure
    adjustedP = correctedP;
    
  case 'Multiple-stage Step-down'
    %this is the step-down version of definition 7 in Benjamini et al. 2006 with L=J
%     q = params.trueNullsEstimationThreshold;
    %k=max{i : for all j<=i there exists l=j so that p(l)<=ql/{m+1?j(1?q)}}.
    %k=max{i : for all j<=i p(j)<=qj/{m+1?j(1?q)}}.
    %this is just the loop version that corresponds to the definition
% % % %     i=1;
% % % %     while  i<numberH0 && all(pData(1:i)<=q*(1:i)'./(numberH0+1-(1:i)'*(1-q)))
% % % %       i = i+1;
% % % %     end
% % % %     k=i;
    %and this a one line version: 
% % % %     k = find(pData>q*(1:numberH0)'./(numberH0+1-(1:numberH0)'*(1-q)),1,'first');
% % % % 
% % % %     correctedP = pData;
% % % %     correctedP(k:end) = 1;
    
    %but I'd rather have an adjustment version:
    %we need to find for each test i, the largest q such that there exists k=max{i : for all j<=i p(j)<=qj/{m+1?j(1?q)}}
    %this is the derivation for q, starting from the previous one-line p-correction:
%     pData>q*(1:numberH0)'./(numberH0+1-(1:numberH0)'*(1-q))
%     pData.*(numberH0+1-(1:numberH0)'*(1-q))>q*(1:numberH0)'
%     pData*(numberH0+1) - pData.*(1:numberH0)'*(1-q) > q*(1:numberH0)'
%     pData*(numberH0+1) - pData.*(1:numberH0)'+ pData.*(1:numberH0)'*q > q*(1:numberH0)'
%     pData*(numberH0+1) - pData.*(1:numberH0)'  > q*(1:numberH0)' - pData.*(1:numberH0)'*q
%     pData*(numberH0+1) - pData.*(1:numberH0)'  > q*((1:numberH0)' - pData.*(1:numberH0)')
%     (pData*(numberH0+1) - pData.*(1:numberH0)')./((1:numberH0)' - pData.*(1:numberH0)')  > q
    adjustedP = (pData*(numberH0+1) - pData.*(1:numberH0)')./((1:numberH0)' - pData.*(1:numberH0)');
    
    %i think we need to enforce monotonicity, but in contrast with the step-up method (see below),
    %it should be done from the smallest to the largest p-value
    for i=2:numberH0
     adjustedP(i) = max(adjustedP(i-1),adjustedP(i));
    end
    adjustedP = min(adjustedP,1); %bound with 1
          
end

adjustedP(sortingIndex) = adjustedP;
adjustedPdata = NaN(sizePdata);
adjustedPdata(isNotNan) = adjustedP;


function qData = linearStepUpFdr(pData,dependenceCorrectionConstant,trueH0Ratio)

numberH0 = length(pData);
%adjustment (does not depend on a pre-chosen threshold)
% it consists in finding, for each p-value, 
% the largest FDR threshold such that the p-value would be considered significant
% it is the inverse operation to p-correction (see commented code below)
% ref: Yekutieli and Benjamini 1999, Journal of Statistical Planning and Inference, 82(1-2) p171
qData = min(pData.*numberH0./(1:numberH0)'.*dependenceCorrectionConstant*trueH0Ratio,1);

%Now, there may be situations where the largest q for a given uncorrected p-value
% is smaller than the largest q for a higher uncorrected p-value 
%(e.g. consecutive identical uncorrected p-values)
% which means that when looking at the map of corrected p-value (=q-values),
% for a given q-value you would reject less hypothesis than in the p-correction method...

%So maybe we should make corrected p-values monotonic in order to have the same number of rejections for the same threshold
%(although i can't find a paper that says that)
%make it monotonic from the end (for each i, q(i) must be smaller than all q between i and n) 
for i=numberH0-1:-1:1
 qData(i) = min(qData(i),qData(i+1));
end