function ea_initrecentpatients(handles,patsub)
% 
%
% USAGE:
%
%    ea_initrecentpatients(handles,patsub)
%
% INPUTS:
%    handles:
%    patsub:
%
% .. AUTHOR:
%       - Andreas Horn, Original file
%       - Ningfei Li, Original file
%       - Daniel Duarte, Documentation

if ~exist('patsub','var')
    patsub='patients';
end
earoot=ea_getearoot;
try
    load([earoot,'common',filesep,'ea_recent',patsub,'.mat']);
catch
    fullrpts={['No recent ',patsub,' found']};
end
save([earoot,'common',filesep,'ea_recent',patsub,'.mat'],'fullrpts');
ea_updaterecentpatients(handles,patsub);