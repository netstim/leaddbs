% enforceMonotonictyResampleFWE.m
%
%        $Id$
%      usage: count = enforceMonotonictyResampleFWE(p,d.indexSortActual,d.indexReorderActual,d.actualIsNotNaN,params)
%         by: julien besle, 
%       date: 14/01/2011
%    purpose: enforces ordered p-value monotonicity in Step-Down bootstrap-based FWE control
%             see algorithm 4.1, step 6, Westfall and Young (1993) p117

function p = enforceMonotonicityResampleFWE(p,d)

if ismember(d.method,{'Step-down','Adaptive Step-down','Adaptive Single-Step'})
  %sort p using provided sort index (sorting of the non-nan actual statistics)
  %note that these values are sorted from the smallest to the largest

  %     reorderedActual = NaN(size(p));                                            % DEBUG
  for i=1:size(p,2)
    sortedP = p(d.indexSortActual{i},i);
  %       sortedActual = actual(d.indexSortActual{i});                                % DEBUG

    %compute consecutive maxima (from largest to smallest p value)
    for j = size(sortedP,1)-1:-1:1
      sortedP(j) = max(sortedP(j),sortedP(j+1));
    end
    p(d.actualIsNotNaN(:,i),i) = sortedP(d.indexReorderActual{i});
%     reorderedActual(d.actualIsNotNaN(:,i),i) = sortedActual(d.indexReorderActual{i});  % DEBUG
  end
  %p(~d.actualIsNotNaN,:)=NaN;
end

