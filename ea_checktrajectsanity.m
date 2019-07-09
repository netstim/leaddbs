function res=checktrajectsanity(trajvector)
% 
%
% USAGE:
%
%    res = checktrajectsanity(trajvector)
%
% INPUT:
%    trajvector:
%
% OUTPUT:
%    res:
%
% .. AUTHOR:
%       - Andreas Horn, Original file
%       - Ning Fey, Original file
%       - Daniel Duarte, Documentation

load('trajvectors');
trajvectors=[trajvectors;trajvector];
save('trajvectors','trajvectors');
res=1;


