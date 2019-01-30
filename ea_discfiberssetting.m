function varargout = ea_discfiberssetting(varargin)
% EA_DISCFIBERSSETTING MATLAB code for ea_discfiberssetting.fig
%      EA_DISCFIBERSSETTING, by itself, creates a new EA_DISCFIBERSSETTING or raises the existing
%      singleton*.
%
%      H = EA_DISCFIBERSSETTING returns the handle to a new EA_DISCFIBERSSETTING or the handle to
%      the existing singleton*.
%
%      EA_DISCFIBERSSETTING('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in EA_DISCFIBERSSETTING.M with the given input arguments.
%
%      EA_DISCFIBERSSETTING('Property','Value',...) creates a new EA_DISCFIBERSSETTING or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ea_discfiberssetting_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ea_discfiberssetting_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ea_discfiberssetting

% Last Modified by GUIDE v2.5 30-Jan-2019 09:32:20

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ea_discfiberssetting_OpeningFcn, ...
                   'gui_OutputFcn',  @ea_discfiberssetting_OutputFcn, ...
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


% --- Executes just before ea_discfiberssetting is made visible.
function ea_discfiberssetting_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ea_discfiberssetting (see VARARGIN)

% Choose default command line output for ea_discfiberssetting
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
    otherwise
        set(handles.showboth, 'Value', 1);
        set(handles.pospredthreshold, 'Enable', 'on');
        set(handles.negpredthreshold, 'Enable', 'on');
end
set(handles.pospredthreshold, 'String', num2str(discfibers.pospredthreshold));
set(handles.negpredthreshold, 'String', num2str(discfibers.negpredthreshold));

% UIWAIT makes ea_discfiberssetting wait for user response (see UIRESUME)
% uiwait(handles.discfiberssetting);


% --- Outputs from this function are returned to the command line.
function varargout = ea_discfiberssetting_OutputFcn(hObject, eventdata, handles) 
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
