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
%       - Ningfei Li, Original file
%       - Daniel Duarte, Documentation

load('trajvectors');
trajvectors=[trajvectors;trajvector];
save('trajvectors','trajvectors');
res=1;


