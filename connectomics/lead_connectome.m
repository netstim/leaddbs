function varargout = lead_connectome(varargin)
% LEADFIGURE MATLAB code for leadfigure.fig
%      LEADFIGURE, by itself, creates a new LEADFIGURE or raises the existing
%      singleton*.
%
%      H = LEADFIGURE returns the handle to a new LEADFIGURE or the handle to
%      the existing singleton*.
%
%      LEADFIGURE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in LEADFIGURE.M with the given input arguments.
%
%      LEADFIGURE('Property','Value',...) creates a new LEADFIGURE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before lead_connectome_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to lead_connectome_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help leadfigure

% Last Modified by GUIDE v2.5 07-Dec-2016 15:43:00

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @lead_connectome_OpeningFcn, ...
                   'gui_OutputFcn',  @lead_connectome_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before leadfigure is made visible.
function lead_connectome_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to leadfigure (see VARARGIN)

earoot=ea_getearoot;
im=imread([earoot,'icons',filesep,'logo_lead_connectome.png']);
image(im);
axes(handles.logoaxes);
axis off;
axis equal;
set(handles.leadfigure,'name','Lead-Connectome','color','w');

% add parcellation atlases to menu:
ll=dir([ea_space([],'labeling'),'*.nii']);
parcellation={};
for lab=1:length(ll)
    [~,n]=fileparts(ll(lab).name);
    parcellation{lab}=n;
end
setappdata(handles.leadfigure,'parcellation',parcellation);
set(handles.parcellation,'String',parcellation);

set(handles.versiontxt,'String',['v',ea_getvsn('local')]);


% add ft methods to menu
cnt=1;
ndir=dir([earoot,'connectomics',filesep,'ea_ft_*.m']);
ftmethod=cell(0);
fdc=cell(0);
for nd=length(ndir):-1:1
    [~,methodf]=fileparts(ndir(nd).name);
    try
        [thisndc,spmvers]=eval([methodf,'(','''prompt''',')']);
        if ismember(spm('ver'),spmvers)
        fdc{cnt}=thisndc;
        ftmethod{cnt}=methodf;
        cnt=cnt+1;
        end
    end
end
setappdata(gcf,'ftmethod',ftmethod);
set(handles.ftmethod,'String',fdc);


% add normmethods to menu
options.prefs=ea_prefs('');
ea_addnormmethods(handles,options,'');

% add recent patients...
ea_initrecentpatients(handles,'subjects');



% update UI:
try
    lc=options.prefs.machine.lc;
catch
    lc=ea_initlcopts(handles.leadfigure);
end
lc2handles(lc,handles);

ea_init_coregmrpopup(handles);

if isempty(varargin) % "standard alone" mode, i.e. not dependend from lead
    isindependent=1;
else
    if strcmp(varargin{1},'dependent')
        isindependent=0;
        ea_makehidelcindep(handles);
    else
        isindependent=1;
    end
end
setappdata(handles.leadfigure,'isindependent',isindependent);

if isindependent
    handles.prod='connectome';
else
    handles.prod='dbs_connectome';
end

if isindependent
    %% add tools menu
    ea_processguiargs(handles,varargin)

    ea_menu_initmenu(handles,{'export','cluster','prefs','transfer','space','methods'});
    
end

ea_firstrun(handles,options);

% Choose default command line output for leadfigure
handles.output = hObject;

% Update handles structure

guidata(hObject, handles);

% UIWAIT makes leadfigure wait for user response (see UIRESUME)

%uiwait(handles.leadfigure);

function ea_makehidelcindep(handles)

set(handles.importpanel,'visible','off');
set(handles.runsavebutn,'String','Save and close');
set(handles.exportcode,'visible','off');

set(handles.openpatientdir,'visible','off');
moveup=165;

pos=get(handles.logoaxes,'Position');
pos(2)=pos(2)-moveup;
set(handles.logoaxes,'Position',pos);

pos=get(handles.titletext,'Position');
pos(2)=pos(2)-moveup;
set(handles.titletext,'Position',pos);

pos=get(handles.versiontxt,'Position');
pos(2)=pos(2)-moveup;
set(handles.versiontxt,'Position',pos);

pos=get(handles.nclinuse,'Position');
pos(2)=pos(2)-moveup;
set(handles.nclinuse,'Position',pos);

pos=get(handles.leadfigure,'Position');
pos(4)=pos(4)-moveup;
set(handles.leadfigure,'Position',pos);


% --- Outputs from this function are returned to the command line.
function varargout = lead_connectome_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1}=hObject;

% --- Executes on selection change in parcellation.
function parcellation_Callback(hObject, eventdata, handles)
% hObject    handle to parcellation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns parcellation contents as cell array
%        contents{get(hObject,'Value')} returns selected item from parcellation


% --- Executes during object creation, after setting all properties.
function parcellation_CreateFcn(hObject, eventdata, handles)
% hObject    handle to parcellation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in compute_CM_func.
function compute_CM_func_Callback(hObject, eventdata, handles)
% hObject    handle to compute_CM_func (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of compute_CM_func


% --- Executes on button press in perf_ft.
function perf_ft_Callback(hObject, eventdata, handles)
% hObject    handle to perf_ft (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of perf_ft


% --- Executes on selection change in ftmethod.
function ftmethod_Callback(hObject, eventdata, handles)
% hObject    handle to ftmethod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns ftmethod contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ftmethod


% --- Executes during object creation, after setting all properties.
function ftmethod_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ftmethod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in compute_CM_struc.
function compute_CM_struc_Callback(hObject, eventdata, handles)
% hObject    handle to compute_CM_struc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of compute_CM_struc


% --- Executes on button press in runsavebutn.
function runsavebutn_Callback(hObject, eventdata, handles)
% hObject    handle to runsavebutn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% check wether this is run or save (independent or dependent mode)

cfig=handles.leadfigure;
options=ea_handles2options(handles);
isindependent=getappdata(handles.leadfigure,'isindependent');
options.uipatdirs=getappdata(cfig,'uipatdir');
ea_savelcopts(handles)

% run execution:

ea_busyaction('on',cfig,'connectome');



if isindependent
ea_run('run',options);
end
ea_busyaction('off',cfig,'connectome');


% --- Executes on button press in degree_centrality.
function degree_centrality_Callback(hObject, eventdata, handles)
% hObject    handle to degree_centrality (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of degree_centrality


% --- Executes on button press in eigenvector_centrality.
function eigenvector_centrality_Callback(hObject, eventdata, handles)
% hObject    handle to eigenvector_centrality (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of eigenvector_centrality


% --- Executes on button press in nodal_efficiency.
function nodal_efficiency_Callback(hObject, eventdata, handles)
% hObject    handle to nodal_efficiency (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of nodal_efficiency


% --- Executes on button press in struc_func_sim.
function struc_func_sim_Callback(hObject, eventdata, handles)
% hObject    handle to struc_func_sim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of struc_func_sim


% --- Executes on button press in compute_GM_func.
function compute_GM_func_Callback(hObject, eventdata, handles)
% hObject    handle to compute_GM_func (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of compute_GM_func


% --- Executes on button press in compute_GM_struc.
function compute_GM_struc_Callback(hObject, eventdata, handles)
% hObject    handle to compute_GM_struc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of compute_GM_struc




function handles=lc2handles(lc,handles)

% General settings

try set(handles.parcellation,'Value',lc.general.parcellationn); end
if get(handles.parcellation,'Value')>length(get(handles.parcellation,'String'))
    set(handles.parcellation,'Value',length(get(handles.parcellation,'String')));
end

% Graph options:
try set(handles.struc_func_sim,'Value',lc.graph.struc_func_sim); end
try set(handles.nodal_efficiency,'Value',lc.graph.nodal_efficiency); end
try set(handles.eigenvector_centrality,'Value',lc.graph.eigenvector_centrality); end
try set(handles.degree_centrality,'Value',lc.graph.degree_centrality); end
try set(handles.fthresh,'String',num2str(lc.graph.fthresh)); end
try set(handles.sthresh,'String',num2str(lc.graph.sthresh)); end


% functional options:
try set(handles.compute_CM_func,'Value',lc.func.compute_CM); end
try set(handles.compute_GM_func,'Value',lc.func.compute_GM); end
try set(handles.TR,'String',num2str(lc.func.prefs.TR)); end


% structural options:
try set(handles.compute_CM_struc,'Value',lc.struc.compute_CM); end

try set(handles.compute_GM_struc,'Value',lc.struc.compute_GM); end
try set(handles.ftmethod,'Value',lc.struc.ft.methodn); end

try set(handles.normalize_fibers,'Value',lc.struc.ft.normalize); end
try set(handles.perf_ft,'Value',lc.struc.ft.do); end


% --- Executes on button press in normalize_fibers.
function normalize_fibers_Callback(hObject, eventdata, handles)
% hObject    handle to normalize_fibers (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of normalize_fibers



function TR_Callback(hObject, eventdata, handles)
% hObject    handle to TR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of TR as text
%        str2double(get(hObject,'String')) returns contents of TR as a double


% --- Executes during object creation, after setting all properties.
function TR_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double


% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function sthresh_Callback(hObject, eventdata, handles)
% hObject    handle to sthresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sthresh as text
%        str2double(get(hObject,'String')) returns contents of sthresh as a double


% --- Executes during object creation, after setting all properties.
function sthresh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sthresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function fthresh_Callback(hObject, eventdata, handles)
% hObject    handle to fthresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fthresh as text
%        str2double(get(hObject,'String')) returns contents of fthresh as a double


% --- Executes during object creation, after setting all properties.
function fthresh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fthresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in normalize_checkbox.
function normalize_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to normalize_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of normalize_checkbox


% --- Executes on selection change in normmethod.
function normmethod_Callback(hObject, eventdata, handles)
% hObject    handle to normmethod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns normmethod contents as cell array
%        contents{get(hObject,'Value')} returns selected item from normmethod
ea_switchnormmethod(handles);


% --- Executes during object creation, after setting all properties.
function normmethod_CreateFcn(hObject, eventdata, handles)
% hObject    handle to normmethod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in normcheck.
function normcheck_Callback(hObject, eventdata, handles)
% hObject    handle to normcheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of normcheck


% --- Executes on button press in dicomcheck.
function dicomcheck_Callback(hObject, eventdata, handles)
% hObject    handle to dicomcheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of dicomcheck
ea_storeui(handles);
ea_deselectall_dicom(handles);

% --- Executes on button press in exportcode.
function exportcode_Callback(hObject, eventdata, handles)
% hObject    handle to exportcode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


ea_savelcopts(handles);


% run execution:

cfig=handles.leadfigure;
ea_busyaction('on',cfig,'connectome');


options=ea_handles2options(handles);
options.macaquemodus=0;
options.uipatdirs=getappdata(handles.leadfigure,'uipatdir');

ea_run('export',options);

ea_busyaction('off',cfig,'connectome');

% --- Executes on button press in patdir_choosebox.
function patdir_choosebox_Callback(hObject, eventdata, handles)
% hObject    handle to patdir_choosebox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ea_busyaction('on',handles.leadfigure,'connectome');
options.prefs=ea_prefs('');
ea_getpatients(options,handles);
ea_busyaction('off',handles.leadfigure,'connectome');

% --- Executes on selection change in recentpts.
function recentpts_Callback(hObject, eventdata, handles)
% hObject    handle to recentpts (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns recentpts contents as cell array
%        contents{get(hObject,'Value')} returns selected item from recentpts
ea_busyaction('on',handles.leadfigure,'connectome');
ea_rcpatientscallback(handles);
ea_busyaction('off',handles.leadfigure,'connectome');


% --- Executes during object creation, after setting all properties.
function recentpts_CreateFcn(hObject, eventdata, handles)
% hObject    handle to recentpts (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in include_lead_connectome_subroutine.
function include_lead_connectome_subroutine_Callback(hObject, eventdata, handles)
% hObject    handle to include_lead_connectome_subroutine (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of include_lead_connectome_subroutine


% --- Executes on button press in coregmr_checkbox.
function coregmr_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to coregmr_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of coregmr_checkbox


% --- Executes on selection change in coregmrpopup.
function coregmrpopup_Callback(hObject, eventdata, handles)
% hObject    handle to coregmrpopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns coregmrpopup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from coregmrpopup


% --- Executes during object creation, after setting all properties.
function coregmrpopup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to coregmrpopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in normsettings.
function normsettings_Callback(hObject, eventdata, handles)
% hObject    handle to normsettings (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
currentNormMethod=getappdata(handles.normsettings,'currentNormMethod');
ea_shownormsettings(currentNormMethod,handles)


% --- Executes on button press in checkfigures.
function checkfigures_Callback(hObject, eventdata, handles)
% hObject    handle to checkfigures (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkfigures


% --- Executes on button press in openpatientdir.
function openpatientdir_Callback(hObject, eventdata, handles)
% hObject    handle to openpatientdir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ea_openpatdir(handles);
