function ea_updaterecentpatients(handles,patsub,nuchosenix)
if ~exist('patsub','var')
    patsub='patients';
end
earoot=ea_getearoot;
load([earoot,'common',filesep,'ea_recentpatients.mat']);
for i=1:length(fullrpts)
    [~,fullrpts{i}]=fileparts(fullrpts{i});
end
fullrpts=[{['Recent ',patsub,':']};fullrpts];
set(handles.recentpts,'String',fullrpts);
if exist('nuchosenix','var')
   set(handles.recentpts,'Value',nuchosenix+1); 
end