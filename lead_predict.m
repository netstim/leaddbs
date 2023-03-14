function varargout = lead_predict(varargin)
% LEAD_PREDICT MATLAB code for lead_predict.fig
%      LEAD_PREDICT, by itself, creates a new LEAD_PREDICT or raises the existing
%      singleton*.
%
%      H = LEAD_PREDICT returns the handle to a new LEAD_PREDICT or the handle to
%      the existing singleton*.
%
%      LEAD_PREDICT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in LEAD_PREDICT.M with the given input arguments.
%
%      LEAD_PREDICT('Property','Value',...) creates a new LEAD_PREDICT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before lead_predict_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to lead_predict_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help lead_predict

% Last Modified by GUIDE v2.5 02-Mar-2023 19:11:49

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @lead_predict_OpeningFcn, ...
                   'gui_OutputFcn',  @lead_predict_OutputFcn, ...
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


% --- Executes just before lead_predict is made visible.
function lead_predict_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   predictmetric line arguments to lead_predict (see VARARGIN)

handles.prod = 'predict';
handles.callingfunction = 'lead_predict';

earoot=ea_getearoot;
im=imread([earoot,'icons',filesep,'logo_lead_predict.png']);
image(im);
axes(handles.logoaxes);
axis off;
axis equal;
set(handles.leadfigure,'name','Lead Predict','color','w');
% homedir=ea_gethome;
%setappdata(handles.leadfigure,'uipatdir',{homedir(1:end-1)});

% add recentpatients patients...
ea_initrecent(handles, 'patients');

ea_processguiargs(handles,varargin)

ea_menu_initmenu(handles,{'cluster','prefs','transfer','vats'},ea_prefs);

[mdl,sf]=ea_genmodlist;
ea_updatemodpopups(mdl,sf,handles)

pmodels=dir([ea_getearoot,'predict',filesep,'ea_predict_*.m']);
for pmod=1:length(pmodels)
    specs=feval(ea_stripext(pmodels(pmod).name),'specs');
    pmods{pmod}=specs.modelname;
    pmodsm{pmod}=ea_stripext(pmodels(pmod).name);
    pmodspecs{pmod}=specs;
end
set(handles.predictionmodel,'String',pmods);
setappdata(handles.predictionmodel,'mfiles',pmodsm);
setappdata(handles.predictionmodel,'specs',pmodspecs);

ea_refreshpredict(handles);

set(handles.versiontxt,'String',['v',ea_getvsn('local')]);

ea_bind_dragndrop(handles.leadfigure, ...
    @(obj,evt) DropFcn(obj,evt,handles), ...
    @(obj,evt) DropFcn(obj,evt,handles));

% Choose default predictmetric line output for lead_predict
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes lead_predict wait for user response (see UIRESUME)
% uiwait(handles.leadfigure);


function ea_refreshpredict(handles)

%modelsel=get(handles.predictionmodel,'String');
%modelsel=modelsel{get(handles.predictionmodel,'Value')};
specs=getappdata(handles.predictionmodel,'specs');
specs=specs{get(handles.predictionmodel,'Value')};

% Coords
set(handles.inccoordinate,'enable',ea_bool2onoff(ismember('Coords',specs.feats)));
set(handles.inccoordinate,'Value',ismember('Coords',specs.feats));
% VTA
set(handles.incvta,'enable',ea_bool2onoff(ismember('VTA',specs.feats)));
set(handles.incvta,'value',ismember('VTA',specs.feats));
% fMRI
set(handles.incfunctional,'enable',ea_bool2onoff(ismember('fMRI',specs.feats)));
set(handles.incfunctional,'value',(ismember('fMRI',specs.feats)));
% dMRI
set(handles.incstructural,'enable',ea_bool2onoff(ismember('dMRI',specs.feats)));
set(handles.incstructural,'value',(ismember('dMRI',specs.feats)));

set(handles.predictmetric,'String',specs.metrics);



% --- Drag and drop callback to load patdir.
function DropFcn(~, event, handles)

switch event.DropType
    case 'file'
        patdir = event.Data;
    case 'string'
        patdir = {event.Data};
end

nonexist = cellfun(@(x) ~exist(x, 'dir'), patdir);
if any(nonexist)
    fprintf('\nExcluded non-existent/invalid folder:\n');
    cellfun(@disp, patdir(nonexist));
    fprintf('\n');
    patdir(nonexist) = [];
end

ea_busyaction('on',handles.leadfigure,'predict');
if ~isempty(patdir)
    ea_load_pts(handles, patdir);
end
ea_busyaction('off',handles.leadfigure,'predict');


% --- Outputs from this function are returned to the predictmetric line.
function varargout = lead_predict_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default predictmetric line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in fiberspopup.
function fiberspopup_Callback(hObject, eventdata, handles)
% hObject    handle to fiberspopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns fiberspopup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from fiberspopup


% --- Executes during object creation, after setting all properties.
function fiberspopup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fiberspopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu5.
function popupmenu5_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu5 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu5


% --- Executes during object creation, after setting all properties.
function popupmenu5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox5.
function checkbox5_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox5


% --- Executes on button press in incstructural.
function incstructural_Callback(hObject, eventdata, handles)
% hObject    handle to incstructural (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of incstructural


% --- Executes on button press in seedbutton.
function seedbutton_Callback(hObject, eventdata, handles)
% hObject    handle to seedbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[seeds,path]=uigetfile({'*.nii','NIfTI';'*.txt','Text';'*.nii.gz','NIfTI'},'Please choose seed definition(s)...','MultiSelect','on');

if iscell(seeds)
    set(hObject,'String',['Multiple (',num2str(length(seeds)),')']);
else
    set(hObject,'String',[seeds]);
    seeds={seeds};
end

for s=1:length(seeds)
   seeds{s}=fullfile(path,seeds{s});
end

setappdata(hObject,'seeds',seeds);


% --- Executes on button press in run_button.
function run_button_Callback(hObject, eventdata, handles)
% hObject    handle to run_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

leadfigure=handles.leadfigure;
ea_busyaction('on',leadfigure,'predict');


options=ea_handles2options(handles);
options.uivatdirs=getappdata(handles.leadfigure,'uipatdir');
options.uipatdirs={''};
options.leadprod = 'predict';

ea_run('run',options);

ea_busyaction('off',leadfigure,'predict');


% --- Executes on button press in exportcode.
function exportcode_Callback(hObject, eventdata, handles)
% hObject    handle to exportcode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

leadfigure=handles.leadfigure;
ea_busyaction('on',leadfigure,'predict');

options=ea_handles2options(handles);
options.uivatdirs=getappdata(handles.leadfigure,'uipatdir');
options.uipatdirs={''};
options.leadprod = 'predict';

ea_run('export',options);

ea_busyaction('off',leadfigure,'predict');


% --- Executes on selection change in predictmetric.
function predictmetric_Callback(hObject, eventdata, handles)
% hObject    handle to predictmetric (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns predictmetric contents as cell array
%        contents{get(hObject,'Value')} returns selected item from predictmetric


% --- Executes during object creation, after setting all properties.
function predictmetric_CreateFcn(hObject, eventdata, handles)
% hObject    handle to predictmetric (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in incfunctional.
function incfunctional_Callback(hObject, eventdata, handles)
% hObject    handle to incfunctional (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of incfunctional


% --- Executes on selection change in fmripopup.
function fmripopup_Callback(hObject, eventdata, handles)
% hObject    handle to fmripopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns fmripopup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from fmripopup


% --- Executes during object creation, after setting all properties.
function fmripopup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fmripopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in odirbutton.
function odirbutton_Callback(hObject, eventdata, handles)
% hObject    handle to odirbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

seeds=getappdata(handles.seedbutton,'seeds');
if ~isempty(seeds) % seeds defined already
   seedbase=fileparts(seeds{1});
else
    seedbase='';
end
odir=uigetdir(seedbase,'Choose output location');
setappdata(hObject,'odir',[odir,filesep]);
set(hObject,'String',odir);


% --- Executes on selection change in strucexportspace.
function strucexportspace_Callback(hObject, eventdata, handles)
% hObject    handle to strucexportspace (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns strucexportspace contents as cell array
%        contents{get(hObject,'Value')} returns selected item from strucexportspace


% --- Executes during object creation, after setting all properties.
function strucexportspace_CreateFcn(hObject, eventdata, handles)
% hObject    handle to strucexportspace (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in omaskbutton.
function omaskbutton_Callback(hObject, eventdata, handles)
% hObject    handle to omaskbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

omask=uigetfile({'*.nii';'*.nii.gz'},'Choose output location');
setappdata(hObject,'omask',[omask]);
set(hObject,'String',omask);


% --- Executes on button press in patdir_choosebox.
function patdir_choosebox_Callback(hObject, eventdata, handles)
% hObject    handle to patdir_choosebox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ea_busyaction('on',handles.leadfigure,'predict');
options.prefs=ea_prefs('');
ea_getpatients(options,handles);
ea_busyaction('off',handles.leadfigure,'predict');

% --- Executes on selection change in recentpatients.
function recentpatients_Callback(hObject, eventdata, handles)
% hObject    handle to recentpatients (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns recentpatients contents as cell array
%        contents{get(hObject,'Value')} returns selected item from recentpatients
ea_busyaction('on',handles.leadfigure,'predict');
ea_recentcallback(handles, 'patients');
ea_busyaction('off',handles.leadfigure,'predict');


% --- Executes during object creation, after setting all properties.
function recentpatients_CreateFcn(hObject, eventdata, handles)
% hObject    handle to recentpatients (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in seeddefpopup.
function seeddefpopup_Callback(hObject, eventdata, handles)
% hObject    handle to seeddefpopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns seeddefpopup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from seeddefpopup

str=get(hObject,'String');
if iscell(str)
    str=str{get(hObject,'Value')};
end

uipatdir=getappdata(handles.leadfigure,'uipatdir');
directory=[uipatdir{1},filesep];
options=ea_handles2options(handles);
options.prefs=ea_prefs;
[mdl,sf]=ea_genmodlist(directory,'nan',options);
ea_updatemodpopups(mdl,sf,handles);


% --- Executes during object creation, after setting all properties.
function seeddefpopup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to seeddefpopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in inccoordinate.
function inccoordinate_Callback(hObject, eventdata, handles)
% hObject    handle to inccoordinate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of inccoordinate


% --- Executes on button press in incvta.
function incvta_Callback(hObject, eventdata, handles)
% hObject    handle to incvta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of incvta


% --- Executes on selection change in predictionmodel.
function predictionmodel_Callback(hObject, eventdata, handles)
% hObject    handle to predictionmodel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns predictionmodel contents as cell array
%        contents{get(hObject,'Value')} returns selected item from predictionmodel
ea_refreshpredict(handles)

% --- Executes during object creation, after setting all properties.
function predictionmodel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to predictionmodel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in usepresentmaps.
function usepresentmaps_Callback(hObject, eventdata, handles)
% hObject    handle to usepresentmaps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of usepresentmaps


% --- Executes on button press in openpatientdir.
function openpatientdir_Callback(hObject, eventdata, handles)
% hObject    handle to openpatientdir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ea_openpatdir(handles);
