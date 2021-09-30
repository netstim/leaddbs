function ea_recentcallback(handles, patsub)

if get(handles.recent,'Value')==1
    return
end

load([ea_getearoot,'common',filesep,'ea_recent',patsub,'.mat']);
if iscell(fullrpts)
    fullrpts=fullrpts(get(handles.recent,'Value')-1);
end

if strcmp(['No recent ' patsub ' found'],fullrpts)
   return
end

if strcmp(patsub,'patients')
    ea_load_pts(handles,fullrpts);
elseif strcmp(patsub,'groups')
    ea_load_group(handles,fullrpts{1,1});
end
