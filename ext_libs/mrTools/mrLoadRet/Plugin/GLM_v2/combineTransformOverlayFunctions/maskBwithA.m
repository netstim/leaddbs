% maskBwithA(A,B,test,newValue)
%
%   masks an overlay B according to a test on the values of overlay A
%     test is an anonymous function that will be applied to A, for example @(x)x<=.05, tests whether A<=.05
%     newValue is the value to replace values in B that DO NOT pass the test (default is nan)
%     
% jb 26/08/2010
%
% $Id$ 
function output = maskBwithA(A,B,test,newValue)

if ~ismember(nargin,[3 4])
   help maskBwithA;
   return
end

if ieNotDefined('newValue')
   newValue = NaN;
end

output = B;

output(~test(A)) = newValue;
