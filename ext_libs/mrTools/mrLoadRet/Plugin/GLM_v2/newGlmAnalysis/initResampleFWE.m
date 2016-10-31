% initResampleFWE.m
%
%        $Id$
%      usage: [indexSortActual,indexReorderActual,d.actualIsNotNaN] = initResampleFWE(actual,params,d.actualIsNotNaN)
%         by: julien besle, 
%       date: 12/01/2011
%    purpose: 
%

function d = initResampleFWE(actual,params,p)

d.method=params.resampleFweMethod;
d.actualIsNotNaN = ~isnan(actual); 
d.numberFalseH0 = zeros(1,size(actual,2));
d.numberTrueH0 = sum(d.actualIsNotNaN,1);

if ismember(d.method,{'Step-down','Adaptive Single-step','Adaptive Step-down'})
  %find the sorting index for actual non-nan T values (independently for each contrast)
  
  actual(~d.actualIsNotNaN) = NaN;
  [sortedActual,d.indexSortActual] = sort(actual,1);
  d.indexSortActual = mat2cell(d.indexSortActual,size(d.indexSortActual,1),ones(size(d.indexSortActual,2),1));
  d.indexReorderActual = cell(size(d.indexSortActual));
  for iTest = 1:length(d.indexSortActual)
    d.indexSortActual{iTest} = d.indexSortActual{iTest}(~isnan(sortedActual(:,iTest))); 
    [temp,d.indexReorderActual{iTest}] = sort(d.indexSortActual{iTest},1);
  end
else
  d.indexSortActual = {};
  d.indexReorderActual = {};
  d.actualIsNotNaN = [];
end

if ismember(d.method,{'Adaptive Single-step','Adaptive Step-down'})
  for iTest = 1:size(p,2)
    numberH0 = length(d.indexSortActual{iTest});
    d.numberTrueH0(iTest) = min(ceil(estimateNumberTrueH0(p(:,iTest),params)),numberH0);
    d.numberFalseH0(iTest) = numberH0 - d.numberTrueH0(iTest);
  end
end



