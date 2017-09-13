function varargout = ea_vatsettings_horn(varargin)
% EA_VATSETTINGS_HORN MATLAB code for ea_vatsettings_horn.fig
%      EA_VATSETTINGS_HORN, by itself, creates a new EA_VATSETTINGS_HORN or raises the existing
%      singleton*.
%
%      H = EA_VATSETTINGS_HORN returns the handle to a new EA_VATSETTINGS_HORN or the handle to
%      the existing singleton*.
%
%      EA_VATSETTINGS_HORN('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in EA_VATSETTINGS_HORN.M with the given input arguments.
%
%      EA_VATSETTINGS_HORN('Property','Value',...) creates a new EA_VATSETTINGS_HORN or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ea_vatsettings_horn_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ea_vatsettings_horn_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ea_vatsettings_horn

% Last Modified by GUIDE v2.5 08-Sep-2017 10:51:02

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ea_vatsettings_horn_OpeningFcn, ...
                   'gui_OutputFcn',  @ea_vatsettings_horn_OutputFcn, ...
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


% --- Executes just before ea_vatsettings_horn is made visible.
function ea_vatsettings_horn_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ea_vatsettings_horn (see VARARGIN)

% Choose default command line output for ea_vatsettings_horn
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

set(handles.setfig,'name','FEM-based VAT model setting');

prefs=ea_prefs('');

set(handles.condgm,'String',num2str(prefs.machine.vatsettings.horn_cgm));
set(handles.condwm,'String',num2str(prefs.machine.vatsettings.horn_cwm));
set(handles.ethresh,'String',num2str(prefs.machine.vatsettings.horn_ethresh));

ea_fillpresetpopups(handles);


function ea_fillpresetpopups(handles)

condv={'Conductivity Presets:',nan,nan
    'General Heuristic (e.g. Buzsaki 2006)',0.33,0.14
    'Frequency Adapted (cf. itis.ethz.ch)',0.0915,0.059};
set(handles.condpresets,'String',condv(:,1));
setappdata(handles.condpresets,'data',condv);


etv={'E-Field Threshold Presets:',nan
    'Heuristic GPI (e.g. Hemm 2005, Vasques 2009)',0.2
    'Heuristic STN (Maedler 2012, Astrom 2014)',0.19
    'Heuristic VIM (Kuncel 2008, Astrom 2014)',0.165
    'D = 2.0 um, 30 us (Astrom 2014)',0.765
    'D = 2.0 um, 60 us (Astrom 2014)',0.457
    'D = 2.0 um, 90 us (Astrom 2014)',0.376
    'D = 2.0 um, 120 us (Astrom 2014)',0.323
    'D = 2.5 um, 30 us (Astrom 2014)',0.504
    'D = 2.5 um, 60 us (Astrom 2014)',0.323
    'D = 2.5 um, 90 us (Astrom 2014)',0.240
    'D = 2.5 um, 120 us (Astrom 2014)',0.310
    'D = 3.0 um, 30 us (Astrom 2014)',0.376
    'D = 3.0 um, 60 us (Astrom 2014)',0.240
    'D = 3.0 um, 90 us (Astrom 2014)',0.185
    'D = 3.0 um, 120 us (Astrom 2014)',0.157
    'D = 3.5 um, 30 us (Astrom 2014)',0.300
    'D = 3.5 um, 60 us (Astrom 2014)',0.185
    'D = 3.5 um, 90 us (Astrom 2014)',0.142
    'D = 3.5 um, 120 us (Astrom 2014)',0.121
    'D = 4.0 um, 30 us (Astrom 2014)',0.240
    'D = 4.0 um, 60 us (Astrom 2014)',0.250
    'D = 4.0 um, 90 us (Astrom 2014)',0.115
    'D = 4.0 um, 120 us (Astrom 2014)',0.096
    'D = 4.5 um, 30 us (Astrom 2014)',0.225
    'D = 4.5 um, 60 us (Astrom 2014)',0.142
    'D = 4.5 um, 90 us (Astrom 2014)',0.107
    'D = 4.5 um, 120 us (Astrom 2014)',0.090
    'D = 5.0 um, 30 us (Astrom 2014)',0.177
    'D = 5.0 um, 60 us (Astrom 2014)',0.110
    'D = 5.0 um, 90 us (Astrom 2014)',0.087
    'D = 5.0 um, 120 us (Astrom 2014)',0.074
    };

set(handles.ethreshpresets,'String',etv(:,1));
setappdata(handles.ethreshpresets,'data',etv);



% UIWAIT makes ea_vatsettings_horn wait for user response (see UIRESUME)


% --- Outputs from this function are returned to the command line.
function varargout = ea_vatsettings_horn_OutputFcn(hObject, eventdata, handles)
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

vatsettings=prefs.machine.vatsettings;
vatsettings.horn_cgm=str2double(get(handles.condgm,'String'));
vatsettings.horn_cwm=str2double(get(handles.condwm,'String'));
vatsettings.horn_ethresh=str2double(get(handles.ethresh,'String'));

ea_setprefs('vatsettings',vatsettings);

delete(handles.setfig);



function ethresh_Callback(hObject, eventdata, handles)
% hObject    handle to ethresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ethresh as text
%        str2double(get(hObject,'String')) returns contents of ethresh as a double


% --- Executes during object creation, after setting all properties.
function ethresh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ethresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in ethreshpresets.
function ethreshpresets_Callback(hObject, eventdata, handles)
% hObject    handle to ethreshpresets (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns ethreshpresets contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ethreshpresets

data=getappdata(hObject,'data');
if isnan(data{get(hObject,'Value'),2})
    return
else
    set(handles.ethresh,'String',num2str(data{get(hObject,'Value'),2}));
end
set(hObject,'value',1);

% --- Executes during object creation, after setting all properties.
function ethreshpresets_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ethreshpresets (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function condgm_Callback(hObject, eventdata, handles)
% hObject    handle to condgm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of condgm as text
%        str2double(get(hObject,'String')) returns contents of condgm as a double


% --- Executes during object creation, after setting all properties.
function condgm_CreateFcn(hObject, eventdata, handles)
% hObject    handle to condgm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in condpresets.
function condpresets_Callback(hObject, eventdata, handles)
% hObject    handle to condpresets (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns condpresets contents as cell array
%        contents{get(hObject,'Value')} returns selected item from condpresets

data=getappdata(hObject,'data');
if isnan(data{get(hObject,'Value'),2})
    return
else
    set(handles.condgm,'String',num2str(data{get(hObject,'Value'),2}));
    set(handles.condwm,'String',num2str(data{get(hObject,'Value'),3}));
end
set(hObject,'value',1);


% --- Executes during object creation, after setting all properties.
function condpresets_CreateFcn(hObject, eventdata, handles)
% hObject    handle to condpresets (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function condwm_Callback(hObject, eventdata, handles)
% hObject    handle to condwm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of condwm as text
%        str2double(get(hObject,'String')) returns contents of condwm as a double


% --- Executes during object creation, after setting all properties.
function condwm_CreateFcn(hObject, eventdata, handles)
% hObject    handle to condwm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
