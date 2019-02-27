function ea_updaterecentpatients(handles,patsub,nuchosenix)
if ~exist('patsub','var')
    patsub='patients';
end
earoot=ea_getearoot;
load([earoot,'common',filesep,'ea_recentpatients.mat']);
for i=1:length(fullrpts)
    [path,foldername]=fileparts(fullrpts{i});
    [~,superfoldername]=fileparts(path);
    fullrpts{i} = [superfoldername filesep foldername];
end
try
    fullrpts=[{['Recent ',patsub,':']};fullrpts];
catch
    fullrpts=[{['Recent ',patsub,':']};fullrpts'];
end
set(handles.recentpts,'String',fullrpts);
if exist('nuchosenix','var')
   set(handles.recentpts,'Value',nuchosenix+1);
end
