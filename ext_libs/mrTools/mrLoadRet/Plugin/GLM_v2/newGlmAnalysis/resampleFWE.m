% resampleFWE.m
%
%        $Id$
%      usage: count = resampleFWE(resampled,actual,d.indexSortActual,params,d)
%         by: julien besle, 
%       date: 12/01/2011
%    purpose: compute Single-Step and Step-Down bootstrap step for bootstrap-based FWE control
%             see algorithm 4.1, step 3, Westfall and Young (1993) p117
%

function count = resampleFWE(resampled,actual,d)

switch(d.method)
  case 'Single-step'
      count = double(repmat(max(resampled,[],1),size(resampled,1), 1)>actual);
    
  case 'Adaptive Single-step'
    count = zeros(size(actual));
    for i=1:size(resampled,2)
      sortedResampled = resampled(d.indexSortActual{i},i);
      maxSortedResampled = sortedResampled;
      for j=1:d.numberFalseH0(i)
        maxSortedResampled(j) = max(sortedResampled(j:j+d.numberTrueH0(i)-1));
      end
      maxSortedResampled(j+1:end) = max(sortedResampled(j+1:end));
      count(d.actualIsNotNaN(:,i),i) = maxSortedResampled(d.indexReorderActual{i});
    end
    count = double(count>actual); %perform the comparison
    count(~d.actualIsNotNaN)=NaN;    
    
  case {'Step-down','Adaptive Step-down'} %Adaptation of step-down (Holm's) procedure
    % In which the probability of the max under complete H0 needs to be estimated on less and less voxels
    %(Holm (1979) Scandinavian Journal of Statistics, Vol. 6, No. 2 , pp. 65-70)
    count = zeros(size(actual));
%     reorderedActual = NaN(size(actual));                                            % DEBUG
    
    for i=1:size(resampled,2)
      sortedResampled = resampled(d.indexSortActual{i},i);
      maxSortedResampled = sortedResampled;
%       sortedActual = actual(d.indexSortActual{i},i);                           % DEBUG

      %compute consecutive maxima such that for increasing statistic values (increasing significance)
      %the max is taken over more and more voxels (with a maximum of d.numberTrueH0 if this has been estimated)
      for j = 2:d.numberTrueH0(i)
        maxSortedResampled(j) = max(sortedResampled(j),maxSortedResampled(j-1));
      end
      %from this point, take only the max over the estimated number of true H0
      for j = d.numberTrueH0(i)+1:length(maxSortedResampled)
        maxSortedResampled(j) = max(sortedResampled(j-d.numberTrueH0(i)+1:j));
      end
      %for each test,compute 1 if the max for this randomization is more than the actual max 
      %temporarily put random maxima in the count variable
      count(d.actualIsNotNaN(:,i),i) = maxSortedResampled(d.indexReorderActual{i});
%       reorderedActual(d.actualIsNotNaN(:,i),i) = sortedActual(d.indexReorderActual{i});  % DEBUG
    end
    count = double(count>actual); %perform the comparison
    count(~d.actualIsNotNaN)=NaN;    
end