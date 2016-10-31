% fdrAdjustedP = fdrAdjustment(p)
%
%   adjusts p values using False Discovery Rate Step-up method
%     
% jb 15/03/2012
%
% $Id: maskAwithB.m 2172 2011-06-20 12:49:44Z julien $ 
function fdrAdjustedP = fdrAdjustment(p)


if ~ismember(nargin,[1])
   help multipleTestsAdjustment;
   return
end

params.fweAdjustment=0;
[~, fdrAdjustedP] = transformStatistic(p,[],params);