% [fdrAdjustedP,fweAdjustedP] = multipleTestsAdjustment(p)
%
%   adjusts p values using False Discovery Rate Step-up method and Hommel Bonferroni correction
%     
% jb 15/03/2012
%
% $Id: maskAwithB.m 2172 2011-06-20 12:49:44Z julien $ 
function [fdrAdjustedP,fweAdjustedP] = multipleTestsAdjustment(p)


if ~ismember(nargin,[1])
   help multipleTestsAdjustment;
   return
end

[~, fdrAdjustedP, fweAdjustedP] = transformStatistic(p);