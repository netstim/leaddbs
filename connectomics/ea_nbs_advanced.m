function varargout = ea_nbs_advanced(varargin)
% EA_NBS_ADVANCED MATLAB code for ea_nbs_advanced.fig
%      EA_NBS_ADVANCED, by itself, creates a new EA_NBS_ADVANCED or raises the existing
%      singleton*.
%
%      H = EA_NBS_ADVANCED returns the handle to a new EA_NBS_ADVANCED or the handle to
%      the existing singleton*.
%
%      EA_NBS_ADVANCED('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in EA_NBS_ADVANCED.M with the given input arguments.
%
%      EA_NBS_ADVANCED('Property','Value',...) creates a new EA_NBS_ADVANCED or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ea_nbs_advanced_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ea_nbs_advanced_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ea_nbs_advanced

% Last Modified by GUIDE v2.5 18-Jun-2016 08:41:22

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ea_nbs_advanced_OpeningFcn, ...
                   'gui_OutputFcn',  @ea_nbs_advanced_OutputFcn, ...
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


% --- Executes just before ea_nbs_advanced is made visible.
function ea_nbs_advanced_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ea_nbs_advanced (see VARARGIN)

% Choose default command line output for ea_nbs_advanced
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ea_nbs_advanced wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = ea_nbs_advanced_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function numpermutations_Callback(hObject, eventdata, handles)
% hObject    handle to numpermutations (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of numpermutations as text
%        str2double(get(hObject,'String')) returns contents of numpermutations as a double


% --- Executes during object creation, after setting all properties.
function numpermutations_CreateFcn(hObject, eventdata, handles)
% hObject    handle to numpermutations (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in component.
function component_Callback(hObject, eventdata, handles)
% hObject    handle to component (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns component contents as cell array
%        contents{get(hObject,'Value')} returns selected item from component


% --- Executes during object creation, after setting all properties.
function component_CreateFcn(hObject, eventdata, handles)
% hObject    handle to component (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in nbsmethod.
function nbsmethod_Callback(hObject, eventdata, handles)
% hObject    handle to nbsmethod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns nbsmethod contents as cell array
%        contents{get(hObject,'Value')} returns selected item from nbsmethod


% --- Executes during object creation, after setting all properties.
function nbsmethod_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nbsmethod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function alpha_Callback(hObject, eventdata, handles)
% hObject    handle to alpha (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of alpha as text
%        str2double(get(hObject,'String')) returns contents of alpha as a double


% --- Executes during object creation, after setting all properties.
function alpha_CreateFcn(hObject, eventdata, handles)
% hObject    handle to alpha (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in savebutton.
function savebutton_Callback(hObject, eventdata, handles)
% hObject    handle to savebutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in exchangebutton.
function exchangebutton_Callback(hObject, eventdata, handles)
% hObject    handle to exchangebutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
