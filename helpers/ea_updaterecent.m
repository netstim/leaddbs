function ea_updaterecent(handles,patsub,nuchosenix)

earoot=ea_getearoot;
load([earoot,'common',filesep,'ea_recent',patsub,'.mat'],'fullrpts');

if strcmp(patsub,'patients')
    for i=1:length(fullrpts)
        [~,fullrpts{i}]=fileparts(fullrpts{i});
    end
end

try
    fullrpts=[{['Recent ',patsub,':']};fullrpts];
catch
    fullrpts=[{['Recent ',patsub,':']};fullrpts'];
end

set(handles.recent,'String',fullrpts);
if exist('nuchosenix','var')
   set(handles.recent,'Value',nuchosenix+1);
end
