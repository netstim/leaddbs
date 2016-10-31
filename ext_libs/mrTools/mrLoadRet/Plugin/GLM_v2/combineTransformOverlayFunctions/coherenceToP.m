% P = coherenceToP(C,n)
%
%   converts Coherence value to p value (under the assumption that time samples are independent)
%         inputs: C = coherence data
%                 n = number of time samples in original time series
%     
% jb 15/03/2012
%
% $Id: maskAwithB.m 2172 2011-06-20 12:49:44Z julien $ 
function P = coherenceToP(C,n)


if ~ismember(nargin,[2])
   help coherenceToP;
   return
end


P = 1-cdf('t',C.*sqrt(n-2)./sqrt(1-C.^2),n-2);