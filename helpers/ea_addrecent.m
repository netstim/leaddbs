function ea_addrecent(handles,uipatdir,patsub,chosenix)

earoot=ea_getearoot;          

load([earoot,'common',filesep,'ea_recent',patsub,'.mat']);
if strcmp(fullrpts,['No recent ',patsub,' found'])
    fullrpts={};
end

if ~exist('chosenix','var')
    try
        chosenix=fullrpts{get(handles.recent,'Value')};
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

ea_updaterecent(handles,patsub,nuchosenix);
