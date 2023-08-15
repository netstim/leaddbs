function varargout = ea_spec2dwrite(varargin)
% EA_SPEC2DWRITE MATLAB code for ea_spec2dwrite.fig
%      EA_SPEC2DWRITE, by itself, creates a new EA_SPEC2DWRITE or raises the existing
%      singleton*.
%
%      H = EA_SPEC2DWRITE returns the handle to a new EA_SPEC2DWRITE or the handle to
%      the existing singleton*.
%
%      EA_SPEC2DWRITE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in EA_SPEC2DWRITE.M with the given input arguments.
%
%      EA_SPEC2DWRITE('Property','Value',...) creates a new EA_SPEC2DWRITE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ea_spec2dwrite_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ea_spec2dwrite_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ea_spec2dwrite

% Last Modified by GUIDE v2.5 18-May-2017 12:54:46

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ea_spec2dwrite_OpeningFcn, ...
                   'gui_OutputFcn',  @ea_spec2dwrite_OutputFcn, ...
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


% --- Executes just before ea_spec2dwrite is made visible.
function ea_spec2dwrite_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ea_spec2dwrite (see VARARGIN)

% Choose default command line output for ea_spec2dwrite
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ea_spec2dwrite wait for user response (see UIRESUME)
% uiwait(handles.ea_spec2dwrite);

% add backdrops

try
    options = varargin{1};
catch
    options.prefs = ea_prefs('');
end

backdrops = ea_assignbackdrop('list',options);

set(handles.tdbackdrop,'String',backdrops);
if get(handles.tdbackdrop,'Value')>length(backdrops)
    set(handles.tdbackdrop,'Value',1);
end

try
    d2=options.prefs.machine.d2;
    ea_options2tdhandles(handles,d2);
end

set(handles.ea_spec2dwrite,'Name','Specify 2D Output options');

ea_ListBoxRenderer(handles.tdbackdrop);


% --- Outputs from this function are returned to the command line.
function varargout = ea_spec2dwrite_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
%varargout{1} = handles.output;


% --- Executes on button press in tdcolorscheck.
function tdcolorscheck_Callback(hObject, eventdata, handles)
% hObject    handle to tdcolorscheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tdcolorscheck


% --- Executes on button press in tdcontourcheck.
function tdcontourcheck_Callback(hObject, eventdata, handles)
% hObject    handle to tdcontourcheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tdcontourcheck



function bbsize_Callback(hObject, eventdata, handles)
% hObject    handle to bbsize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of bbsize as text
%        str2double(get(hObject,'String')) returns contents of bbsize as a double


% --- Executes during object creation, after setting all properties.
function bbsize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bbsize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in tdlabelcheck.
function tdlabelcheck_Callback(hObject, eventdata, handles)
% hObject    handle to tdlabelcheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tdlabelcheck


% --- Executes on button press in tdcontourcolor.
function tdcontourcolor_Callback(hObject, eventdata, handles)
% hObject    handle to tdcontourcolor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
tdcol = ea_uisetcolor;
setappdata(hObject,'color',tdcol);

% --- Executes on button press in savebutton.
function savebutton_Callback(hObject, eventdata, handles)
% hObject    handle to savebutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
d2=ea_tdhandles2options(handles);
ea_setprefs('d2',d2);
delete(handles.ea_spec2dwrite);


% --- Executes on selection change in tdbackdrop.
function tdbackdrop_Callback(hObject, eventdata, handles)
% hObject    handle to tdbackdrop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns tdbackdrop contents as cell array
%        contents{get(hObject,'Value')} returns selected item from tdbackdrop


% --- Executes during object creation, after setting all properties.
function tdbackdrop_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tdbackdrop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in tdfidcheck.
function tdfidcheck_Callback(hObject, eventdata, handles)
% hObject    handle to tdfidcheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tdfidcheck
