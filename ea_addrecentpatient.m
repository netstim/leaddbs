function ea_addrecentpatient(handles,uipatdir,patsub,chosenix)
earoot=ea_getearoot;
if ~exist('patsub','var')
    patsub='patients';
end
load([earoot,'common',filesep,'ea_recentpatients.mat']);
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
fullrpts=[uipatdir;fullrpts];    
end

[fullrpts]=unique(fullrpts,'stable');
if length(fullrpts)>10
    
   fullrpts=fullrpts(1:10);
end
[~,nuchosenix]=ismember(chosenix,fullrpts);
save([earoot,'common',filesep,'ea_recentpatients.mat'],'fullrpts');

ea_updaterecentpatients(handles,patsub,nuchosenix);
