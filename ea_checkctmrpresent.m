function modality=ea_checkctmrpresent(varargin)

if ischar(varargin{1})
    uipatdir{1}=varargin{1};
else % call from lead directly:
    handles=varargin{1};
    uipatdir=getappdata(handles.leadfigure,'uipatdir');
end

% check for MR-files
[~,patientname]=fileparts(uipatdir{1});
prefs=ea_prefs(patientname);

modality=zeros(2,1);
% check for unnormalized images first:
% CT
if exist([uipatdir{1},filesep,prefs.rawctnii_unnormalized],'file')
try    set(handles.statusone,'String','Raw CT file found. Please coregister CT to preoperative MR first.'); end
    modality(2)=1;
end
if exist([uipatdir{1},filesep,prefs.ctnii_coregistered],'file')
try    set(handles.statusone,'String','Coregistered CT file found. Set normalize option.'); end
    modality(2)=1;
end
% MR
if exist([uipatdir{1},filesep,prefs.tranii_unnormalized],'file')
try    set(handles.statusone,'String','Unnormalized MR-volumes found. Set normalize option.'); end
    modality(1)=1;
end
% MR & CT
if exist([uipatdir{1},filesep,prefs.tranii_unnormalized],'file') && exist([uipatdir{1},filesep,prefs.cornii_unnormalized],'file') && exist([uipatdir{1},filesep,prefs.rawctnii_unnormalized],'file')
try    set(handles.statusone,'String','Unnormalized MR-volumes and raw CT volume found. Set normalize option or coregister CT.'); end
    modality(1:2)=1;
end
if exist([uipatdir{1},filesep,prefs.tranii_unnormalized],'file') && exist([uipatdir{1},filesep,prefs.cornii_unnormalized],'file') && exist([uipatdir{1},filesep,prefs.ctnii_coregistered],'file')
try    set(handles.statusone,'String','Unnormalized MR-volumes and coregistered CT file found. Set normalize option.'); end
    modality(1:2)=1;
end


% Now check for normalized images:
if exist([uipatdir{1},filesep,prefs.ctnii],'file')
    try
    load([uipatdir{1},filesep,'ea_normmethod_applied.mat'])
    nmused=eval([norm_method_applied,'(''check'')']);
try    set(handles.statusone,'String',['CT-volumes have been normalized with: ',nmused]); end
    catch
try    set(handles.statusone,'String','Normalized CT-volumes found.'); end
    end
    modality(2)=1;
end
if exist([uipatdir{1},filesep,prefs.tranii],'file') && exist([uipatdir{1},filesep,prefs.cornii],'file')
    try
    load([uipatdir{1},filesep,'ea_normmethod_applied.mat'])
    nmused=eval([norm_method_applied,'(''check'')']);
try    set(handles.statusone,'String',['MR-volumes have been normalized with: ',nmused]); end
    catch
try    set(handles.statusone,'String','Normalized MR-volumes found.'); end
    end
    modality(1)=1;
end
if exist([uipatdir{1},filesep,prefs.tranii],'file') && exist([uipatdir{1},filesep,prefs.cornii],'file') && exist([uipatdir{1},filesep,prefs.ctnii],'file')
 try   set(handles.statusone,'String','Normalized MR- and CT-volumes found.'); end
    modality(1:2)=1;
end