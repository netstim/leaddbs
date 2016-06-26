function varargout = lead_connectome(varargin)
% LEAD_CONNECTOME MATLAB code for lead_connectome.fig
%      LEAD_CONNECTOME, by itself, creates a new LEAD_CONNECTOME or raises the existing
%      singleton*.
%
%      H = LEAD_CONNECTOME returns the handle to a new LEAD_CONNECTOME or the handle to
%      the existing singleton*.
%
%      LEAD_CONNECTOME('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in LEAD_CONNECTOME.M with the given input arguments.
%
%      LEAD_CONNECTOME('Property','Value',...) creates a new LEAD_CONNECTOME or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before lead_connectome_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to lead_connectome_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help lead_connectome

% Last Modified by GUIDE v2.5 14-May-2015 19:47:03

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


% --- Executes just before lead_connectome is made visible.
function lead_connectome_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to lead_connectome (see VARARGIN)

earoot=[ea_getearoot];
im=imread('ea_logo.png');
image(im);
axis off;
axis equal;
set(gcf,'name','Welcome to LEAD-DBS','color','w');

% add parcellation atlases to menu:
ll=dir([ea_getearoot,'templates',filesep,'labeling',filesep,'*.nii']);
for lab=1:length(ll)
    [~,n]=fileparts(ll(lab).name);
    parcellation{lab}=n;
end
setappdata(gcf,'parcellation',parcellation);
set(handles.parcellation,'String',parcellation);


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



% update UI:
try
    lc=load([ea_getearoot,'connectomics',filesep,'lc_options.mat']);
catch
    lc=ea_initlcopts(handles.lead_connectome);
end
lc2handles(lc,handles);


% Choose default command line output for lead_connectome
handles.output = hObject;

% Update handles structure

guidata(hObject, handles);

% UIWAIT makes lead_connectome wait for user response (see UIRESUME)

%uiwait(handles.lead_connectome);


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


% --- Executes on button press in savebutn.
function savebutn_Callback(hObject, eventdata, handles)
% hObject    handle to savebutn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

lc_options=handles2lc(handles);
save([ea_getearoot,'connectomics',filesep,'lc_options.mat'],'-struct','lc_options');
delete(handles.lead_connectome);
uiresume


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


function lc=handles2lc(handles)


% General settings

lc.general.parcellation=getappdata(gcf,'parcellation');
lc.general.parcellation=lc.general.parcellation{get(handles.parcellation,'Value')};
lc.general.parcellationn=get(handles.parcellation,'Value');


% Graph options:
lc.graph.struc_func_sim=get(handles.struc_func_sim,'Value');
lc.graph.nodal_efficiency=get(handles.nodal_efficiency,'Value');
lc.graph.eigenvector_centrality=get(handles.eigenvector_centrality,'Value');
lc.graph.degree_centrality=get(handles.degree_centrality,'Value');
lc.graph.fthresh=str2double(get(handles.fthresh,'String'));
lc.graph.sthresh=str2double(get(handles.sthresh,'String'));


% functional options:
lc.func.compute_CM=get(handles.compute_CM_func,'Value');
lc.func.compute_GM=get(handles.compute_GM_func,'Value');
lc.func.prefs.TR=str2double(get(handles.TR,'String'));


% structural options:
lc.struc.compute_CM=get(handles.compute_CM_struc,'Value');
lc.struc.compute_GM=get(handles.compute_GM_struc,'Value');
lc.struc.ft.method=getappdata(gcf,'ftmethod');
lc.struc.ft.method=lc.struc.ft.method{get(handles.ftmethod,'Value')};
lc.struc.ft.methodn=get(handles.ftmethod,'Value');
lc.struc.ft.do=get(handles.perf_ft,'Value');
lc.struc.ft.normalize=get(handles.normalize_fibers,'Value');











function handles=lc2handles(lc,handles)

% General settings

try set(handles.parcellation,'Value',lc.general.parcellationn); end


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




