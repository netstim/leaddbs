function varargout = ea_discfibers_setting(varargin)
% EA_DISCFIBERS_SETTING MATLAB code for ea_discfibers_setting.fig
%      EA_DISCFIBERS_SETTING, by itself, creates a new EA_DISCFIBERS_SETTING or raises the existing
%      singleton*.
%
%      H = EA_DISCFIBERS_SETTING returns the handle to a new EA_DISCFIBERS_SETTING or the handle to
%      the existing singleton*.
%
%      EA_DISCFIBERS_SETTING('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in EA_DISCFIBERS_SETTING.M with the given input arguments.
%
%      EA_DISCFIBERS_SETTING('Property','Value',...) creates a new EA_DISCFIBERS_SETTING or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ea_discfibers_setting_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ea_discfibers_setting_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ea_discfibers_setting

% Last Modified by GUIDE v2.5 17-Feb-2019 15:46:38

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ea_discfibers_setting_OpeningFcn, ...
                   'gui_OutputFcn',  @ea_discfibers_setting_OutputFcn, ...
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


% --- Executes just before ea_discfibers_setting is made visible.
function ea_discfibers_setting_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ea_discfibers_setting (see VARARGIN)

% Choose default command line output for ea_discfibers_setting
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% select last selection
prefs = ea_prefs('');
discfibers = prefs.machine.lg.discfibers;
set(handles.connthreshold, 'String', num2str(discfibers.connthreshold));
switch discfibers.showfibersset
    case 'positive'
        set(handles.showposonly, 'Value', 1);
        set(handles.pospredthreshold, 'Enable', 'on');
        set(handles.negpredthreshold, 'Enable', 'off');
    case 'negative'
        set(handles.shownegonly, 'Value', 1);
        set(handles.pospredthreshold, 'Enable', 'off');
        set(handles.negpredthreshold, 'Enable', 'on');
    case 'both'
        set(handles.showboth, 'Value', 1);
        set(handles.pospredthreshold, 'Enable', 'on');
        set(handles.negpredthreshold, 'Enable', 'on');
end
set(handles.pospredthreshold, 'String', num2str(discfibers.pospredthreshold));
set(handles.negpredthreshold, 'String', num2str(discfibers.negpredthreshold));
set(handles.statmetric,'Value',discfibers.statmetric);


if get(handles.statmetric,'Value')>1
    set(handles.connthreshold,'enable','off');
else
        set(handles.connthreshold,'enable','on');
end

% UIWAIT makes ea_discfibers_setting wait for user response (see UIRESUME)
% uiwait(handles.discfiberssetting);


% --- Outputs from this function are returned to the command line.
function varargout = ea_discfibers_setting_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


function connthreshold_Callback(hObject, eventdata, handles)
% hObject    handle to connthreshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of connthreshold as text
%        str2double(get(hObject,'String')) returns contents of connthreshold as a double


% --- Executes during object creation, after setting all properties.
function connthreshold_CreateFcn(hObject, eventdata, handles)
% hObject    handle to connthreshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function pospredthreshold_Callback(hObject, eventdata, handles)
% hObject    handle to pospredthreshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pospredthreshold as text
%        str2double(get(hObject,'String')) returns contents of pospredthreshold as a double


% --- Executes during object creation, after setting all properties.
function pospredthreshold_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pospredthreshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in savesetting.
function savesetting_Callback(hObject, eventdata, handles)
% hObject    handle to savesetting (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

prefs=ea_prefs('');
discfibers = prefs.machine.lg.discfibers;
discfibers.connthreshold = str2double(get(handles.connthreshold,'String'));
switch get(get(handles.showfiberssetpanel, 'SelectedObject'), 'Tag')
    case 'showposonly'
        discfibers.showfibersset = 'positive';
    case 'shownegonly'
        discfibers.showfibersset = 'negative';
    case 'showboth'
        discfibers.showfibersset = 'both';
end
discfibers.pospredthreshold = str2double(get(handles.pospredthreshold,'String'));
discfibers.negpredthreshold = str2double(get(handles.negpredthreshold,'String'));
discfibers.statmetric=get(handles.statmetric,'Value');
ea_setprefs('lg.discfibers', discfibers);

delete(handles.discfiberssetting);


% --- Executes on button press in showposonly.
function showposonly_Callback(hObject, eventdata, handles)
% hObject    handle to showposonly (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of showposonly
set(handles.pospredthreshold, 'Enable', 'on');
set(handles.negpredthreshold, 'Enable', 'off');


% --- Executes on button press in shownegonly.
function shownegonly_Callback(hObject, eventdata, handles)
% hObject    handle to shownegonly (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of shownegonly
set(handles.pospredthreshold, 'Enable', 'off');
set(handles.negpredthreshold, 'Enable', 'on');


% --- Executes on button press in showboth.
function showboth_Callback(hObject, eventdata, handles)
% hObject    handle to showboth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of showboth
set(handles.pospredthreshold, 'Enable', 'on');
set(handles.negpredthreshold, 'Enable', 'on');



function negpredthreshold_Callback(hObject, eventdata, handles)
% hObject    handle to negpredthreshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of negpredthreshold as text
%        str2double(get(hObject,'String')) returns contents of negpredthreshold as a double


% --- Executes during object creation, after setting all properties.
function negpredthreshold_CreateFcn(hObject, eventdata, handles)
% hObject    handle to negpredthreshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in statmetric.
function statmetric_Callback(hObject, eventdata, handles)
% hObject    handle to statmetric (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns statmetric contents as cell array
%        contents{get(hObject,'Value')} returns selected item from statmetric
if get(handles.statmetric,'Value')>1
    set(handles.connthreshold,'enable','off');
else
    set(handles.connthreshold,'enable','on');
end

% --- Executes during object creation, after setting all properties.
function statmetric_CreateFcn(hObject, eventdata, handles)
% hObject    handle to statmetric (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
