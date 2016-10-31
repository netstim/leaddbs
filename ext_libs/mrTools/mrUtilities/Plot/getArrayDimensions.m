
function [i,j]=getArrayDimensions(n,maxRatio)

%
%       [i,j]=getArrayDimensions(n,maxRatio)
%
%        $Id: getArrayDimensions.m 2134 2011-05-31 08:49:17Z julien $
%             jb 01/06/2011
%
%           returns a number of rows and columns for a matrix of n elements 
%           with a ratio rows/columns approaching but less than maxRatio

if ieNotDefined('maxRatio')
  maxRatio = .75;
end

for i=1:n
  if i/ceil(n/i)>maxRatio
    i = max(1,i-1);
    j = ceil(n/i);
    return
  end
end

j=1;

end
