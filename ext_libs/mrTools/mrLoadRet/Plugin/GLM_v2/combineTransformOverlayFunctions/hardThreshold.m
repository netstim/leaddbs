% hardThreshold(B,A,test,newValue)
%
%   masks an overlay A according to a test on its values
%     test is an anonymous function that will be applied to B, for example @(x)x<=.05, tests whether A<=.05
%     newValue is the value to replace values in A that DO NOT pass the test (default is nan)
%     
% jb 27/08/2011
%
% $Id: maskAwithB.m 1950 2010-12-18 10:12:48Z julien $ 
function output = hardThreshold(A,test,newValue)


if ~ismember(nargin,[2 3])
   help maskAwithB;
   return
end

if ieNotDefined('newValue')
   newValue = NaN;
end

output = maskBwithA(A,A,test,newValue);