% Z = PtoZ(P)
%
%   converts P value to Z value 
%     
% jb 15/03/2012
%
% $Id: maskAwithB.m 2172 2011-06-20 12:49:44Z julien $ 
function Z = PtoZ(P)


if ~ismember(nargin,[1])
   help coherenceToP;
   return
end

params.fweAdjustment= 0;
params.fdrAdjustment= 0;
params.testOutput = 'Z value';
Z = transformStatistic(P,[],params);