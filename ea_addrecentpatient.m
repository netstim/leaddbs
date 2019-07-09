function ea_addrecentpatient(handles,uipatdir,patsub,chosenix)
% 
%
% USAGE:
%
%    ea_addrecentpatient(handles,uipatdir,patsub,chosenix)
%
% INPUT:
%    handles:
%    uipatdir:
%    patsub:
%    chosenix:
%
% .. AUTHOR:
%       - Andreas Horn, Original file
%       - Ning Fey, Original file
%       - Daniel Duarte, Documentation

earoot=ea_getearoot;          

load([earoot,'common',filesep,'ea_recent',patsub,'.mat']);
if strcmp(fullrpts,['No recent ',patsub,' found'])
    fullrpts={};
end

if ~exist('chosenix','var')
    try
        chosenix=fullrpts{get(handles.recentpts,'Value')};
    catch
        chosenix=['Recent ',patsub,':'];
    end
end

try
    fullrpts=[uipatdir';fullrpts];
catch % calls from lead_group could end up transposed
    try
        fullrpts=[uipatdir;fullrpts];
    catch
        fullrpts=[uipatdir;fullrpts'];
    end
end

[fullrpts]=unique(fullrpts,'stable');
if length(fullrpts)>10

   fullrpts=fullrpts(1:10);
end
[~,nuchosenix]=ismember(chosenix,fullrpts);
save([earoot,'common',filesep,'ea_recent',patsub,'.mat'],'fullrpts');

ea_updaterecentpatients(handles,patsub,nuchosenix);
