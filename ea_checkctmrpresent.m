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
    status = 'Raw CT file found. Please coregister CT to preoperative MR first.';
    modality(2)=1;
end
if exist([uipatdir{1},filesep,prefs.ctnii_coregistered],'file')
    status = 'Coregistered CT file found. Set Normalization option.';
    modality(2)=1;
end

% MR
if exist([uipatdir{1},filesep,prefs.tranii_unnormalized],'file')
    status = 'Unnormalized MR-volumes found. Set Normalization option.';
    modality(1)=1;
end

% MR & CT
if exist([uipatdir{1},filesep,prefs.tranii_unnormalized],'file') && exist([uipatdir{1},filesep,prefs.rawctnii_unnormalized],'file')
    status = 'Unnormalized MR-volumes and raw CT volume found. Set Normalization option or Coregister CT.';
    modality(1:2)=1;
end
if exist([uipatdir{1},filesep,prefs.tranii_unnormalized],'file') && exist([uipatdir{1},filesep,prefs.ctnii_coregistered],'file')
    status = 'Unnormalized MR-volumes and coregistered CT file found. Set Normalization option.';
    modality(1:2)=1;
end

% Now check for normalized images:
if exist([uipatdir{1},filesep,prefs.ctnii],'file')
    try
        load([uipatdir{1},filesep,'ea_normmethod_applied.mat'])
        nmused=eval([norm_method_applied{1},'(''check'')']);
        status = ['CT-volumes have been normalized with: ',nmused];
    catch
        status = 'Normalized CT-volumes found.';
    end
    modality(2)=1;
end

if exist([uipatdir{1},filesep,prefs.tranii],'file')
    try
        load([uipatdir{1},filesep,'ea_normmethod_applied.mat'])
        nmused=eval([norm_method_applied{1},'(''check'')']);
        status = ['MR-volumes have been normalized with: ',nmused];
    catch
        status = 'Normalized MR-volumes found.';
    end
    modality(1)=1;
end

if exist([uipatdir{1},filesep,prefs.tranii],'file') && exist([uipatdir{1},filesep,prefs.ctnii],'file')
    status = 'Normalized MR- and CT-volumes found.';
    modality(1:2)=1;
end
try % not available when calling from group
set(handles.statusone,'String',status);
set(handles.statusone,'TooltipString',status);
end