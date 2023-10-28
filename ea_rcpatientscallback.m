function ea_rcpatientscallback(handles, patsub)

if get(handles.recentpts,'Value')==1
    return
end

load([ea_prefsdir, filesep, 'ea_recent', patsub, '.mat']);
if iscell(fullrpts)
    fullrpts=fullrpts(get(handles.recentpts,'Value')-1);
end

if strcmp(['No recent ' patsub ' found'],fullrpts)
   return
end

if strcmp(patsub,'patients')
    ea_load_pts(handles,fullrpts);
elseif strcmp(patsub,'groups')
    ea_load_group(handles,fullrpts{1,1});
end
