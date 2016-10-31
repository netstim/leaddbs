% maskAwithB(A,B,test,newValue)
%
%   masks an overlay A according to a test on the values of overlay B
%     test is an anonymous function that will be applied to B, for example @(x)x<=.05, tests whether B<=.05
%     newValue is the value to replace values in A that DO NOT pass the test (default is nan)
%     
% jb 27/08/2010
%
% $Id$ 
function output = maskAwithB(A,B,test,newValue)


if ~ismember(nargin,[3 4])
   help maskAwithB;
   return
end

if ieNotDefined('newValue')
   newValue = NaN;
end

output = maskBwithA(B,A,test,newValue);