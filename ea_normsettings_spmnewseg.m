function varargout = ea_normsettings_spmnewseg(varargin)
% EA_NORMSETTINGS_SPMNEWSEG MATLAB code for ea_normsettings_spmnewseg.fig
%      EA_NORMSETTINGS_SPMNEWSEG, by itself, creates a new EA_NORMSETTINGS_SPMNEWSEG or raises the existing
%      singleton*.
%
%      H = EA_NORMSETTINGS_SPMNEWSEG returns the handle to a new EA_NORMSETTINGS_SPMNEWSEG or the handle to
%      the existing singleton*.
%
%      EA_NORMSETTINGS_SPMNEWSEG('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in EA_NORMSETTINGS_SPMNEWSEG.M with the given input arguments.
%
%      EA_NORMSETTINGS_SPMNEWSEG('Property','Value',...) creates a new EA_NORMSETTINGS_SPMNEWSEG or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ea_normsettings_spmnewseg_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ea_normsettings_spmnewseg_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ea_normsettings_spmnewseg

% Last Modified by GUIDE v2.5 28-Nov-2018 17:24:35

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ea_normsettings_spmnewseg_OpeningFcn, ...
                   'gui_OutputFcn',  @ea_normsettings_spmnewseg_OutputFcn, ...
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


% --- Executes just before ea_normsettings_spmnewseg is made visible.
function ea_normsettings_spmnewseg_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ea_normsettings_spmnewseg (see VARARGIN)

% Choose default command line output for ea_normsettings_spmnewseg
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

set(handles.setfig,'Name','Normalization Settings');
set(handles.titletext,'String','SPM Segment Defaults');




% list presets

earoot=ea_getearoot;

% select last selection
prefs=ea_prefs('');
set(handles.scalereg,'String',num2str(prefs.machine.normsettings.spmnewseg_scalereg));

% UIWAIT makes ea_normsettings_spmnewseg wait for user response (see UIRESUME)


% --- Outputs from this function are returned to the command line.
function varargout = ea_normsettings_spmnewseg_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure

% --- Executes on selection change in pcpopup.
function pcpopup_Callback(hObject, eventdata, handles)
% hObject    handle to pcpopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns pcpopup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pcpopup
selpc=get(hObject,'String');
selpc=selpc{get(hObject,'Value')};
if strcmp(selpc,'User-Defined')
   uipeerdir=ea_getpatients;
   setappdata(hObject,'peersetcell',uipeerdir);
end


% --- Executes during object creation, after setting all properties.
function pcpopup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pcpopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in savebutn.
function savebutn_Callback(hObject, eventdata, handles)
% hObject    handle to savebutn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
prefs=ea_prefs('');
normsettings=prefs.machine.normsettings;
normsettings.spmnewseg_scalereg=str2double(get(handles.scalereg,'String'));

ea_setprefs('normsettings',normsettings);

delete(handles.setfig);


% --- Executes on button press in scrf.
function scrf_Callback(hObject, eventdata, handles)
% hObject    handle to scrf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of scrf


% --- Executes on selection change in metric.
function metric_Callback(hObject, eventdata, handles)
% hObject    handle to metric (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns metric contents as cell array
%        contents{get(hObject,'Value')} returns selected item from metric


% --- Executes during object creation, after setting all properties.
function metric_CreateFcn(hObject, eventdata, handles)
% hObject    handle to metric (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in strategy.
function strategy_Callback(hObject, eventdata, handles)
% hObject    handle to strategy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns strategy contents as cell array
%        contents{get(hObject,'Value')} returns selected item from strategy


% --- Executes during object creation, after setting all properties.
function strategy_CreateFcn(hObject, eventdata, handles)
% hObject    handle to strategy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function numcores_Callback(hObject, eventdata, handles)
% hObject    handle to numcores (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of numcores as text
%        str2double(get(hObject,'String')) returns contents of numcores as a double


% --- Executes during object creation, after setting all properties.
function numcores_CreateFcn(hObject, eventdata, handles)
% hObject    handle to numcores (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in skullstrip.
function skullstrip_Callback(hObject, eventdata, handles)
% hObject    handle to skullstrip (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of skullstrip



function scalereg_Callback(hObject, eventdata, handles)
% hObject    handle to scalereg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of scalereg as text
%        str2double(get(hObject,'String')) returns contents of scalereg as a double


% --- Executes during object creation, after setting all properties.
function scalereg_CreateFcn(hObject, eventdata, handles)
% hObject    handle to scalereg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
