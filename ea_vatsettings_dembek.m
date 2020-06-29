function varargout = ea_vatsettings_dembek(varargin)
% EA_VATSETTINGS_DEMBEK MATLAB code for ea_vatsettings_dembek.fig
%      EA_VATSETTINGS_DEMBEK, by itself, creates a new EA_VATSETTINGS_DEMBEK or raises the existing
%      singleton*.
%
%      H = EA_VATSETTINGS_DEMBEK returns the handle to a new EA_VATSETTINGS_DEMBEK or the handle to
%      the existing singleton*.
%
%      EA_VATSETTINGS_DEMBEK('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in EA_VATSETTINGS_DEMBEK.M with the given input arguments.
%
%      EA_VATSETTINGS_DEMBEK('Property','Value',...) creates a new EA_VATSETTINGS_DEMBEK or raises the
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

% Last Modified by GUIDE v2.5 04-May-2018 12:25:41

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ea_vatsettings_dembek_OpeningFcn, ...
                   'gui_OutputFcn',  @ea_vatsettings_dembek_OutputFcn, ...
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
function ea_vatsettings_dembek_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ea_vatsettings_dembek (see VARARGIN)

% Choose default command line output for ea_vatsettings_dembek
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

set(handles.setfig,'name','FEM-based VAT model setting');

prefs=ea_prefs('');

% set(handles.condgm,'String',num2str(prefs.machine.vatsettings.horn_cgm));
% set(handles.condwm,'String',num2str(prefs.machine.vatsettings.horn_cwm));
set(handles.ethreshtext,'String',num2str(prefs.machine.vatsettings.dembek_ethresh));  
set(handles.ethreshtext,'Value',prefs.machine.vatsettings.dembek_ethreshpw);  
set(handles.pwtext,'String',num2str(prefs.machine.vatsettings.dembek_pw));  

ea_fillpresetpopups(handles);


function ea_fillpresetpopups(handles)

% condv={'Conductivity Presets:',nan,nan
%     'General Heuristic (e.g. Buzsaki 2006)',0.33,0.14
%     'Frequency Adapted (cf. itis.ethz.ch)',0.0915,0.059};
% set(handles.pwpresets,'String',condv(:,1));
% setappdata(handles.pwpresets,'data',condv);


pwp={'Stimulation Pulse Width:',nan
    '30 us',30
    '60 us',60
    '90 us',90
    '120 us',120
    };
set(handles.pwpresets,'String',pwp(:,1));
setappdata(handles.pwpresets,'data',pwp);

etv={'E-Field Threshold Presets:',nan,nan
    'General Heuristic (e.g. Hemm 2005. Vasques 2009, Astrom 2009, Horn 2017)',0.2,60
    'Heuristic GPI (e.g. Hemm 2005, Vasques 2009)',0.2,60
    'Heuristic STN (Maedler 2012, Astrom 2014)',0.19,60
    'Heuristic VIM (Kuncel 2008, Astrom 2014)',0.165,90
    'D = 2.0 um, 30 us (Astrom 2014)',0.765,30
    'D = 2.0 um, 60 us (Astrom 2014)',0.457,60
    'D = 2.0 um, 90 us (Astrom 2014)',0.376,90
    'D = 2.0 um, 120 us (Astrom 2014)',0.323,120
    'D = 2.5 um, 30 us (Astrom 2014)',0.504,30
    'D = 2.5 um, 60 us (Astrom 2014)',0.323,60
    'D = 2.5 um, 90 us (Astrom 2014)',0.240,90
    'D = 2.5 um, 120 us (Astrom 2014)',0.310,120
    'D = 3.0 um, 30 us (Astrom 2014)',0.376,30
    'D = 3.0 um, 60 us (Astrom 2014)',0.240,60
    'D = 3.0 um, 90 us (Astrom 2014)',0.185,90
    'D = 3.0 um, 120 us (Astrom 2014)',0.157,120
    'D = 3.5 um, 30 us (Astrom 2014)',0.300,30
    'D = 3.5 um, 60 us (Astrom 2014)',0.185,60
    'D = 3.5 um, 90 us (Astrom 2014)',0.142,90
    'D = 3.5 um, 120 us (Astrom 2014)',0.121,120
    'D = 4.0 um, 30 us (Astrom 2014)',0.240,30
    'D = 4.0 um, 60 us (Astrom 2014)',0.150,60
    'D = 4.0 um, 90 us (Astrom 2014)',0.115,90
    'D = 4.0 um, 120 us (Astrom 2014)',0.096,120
    'D = 4.5 um, 30 us (Astrom 2014)',0.225,30
    'D = 4.5 um, 60 us (Astrom 2014)',0.142,60
    'D = 4.5 um, 90 us (Astrom 2014)',0.107,90
    'D = 4.5 um, 120 us (Astrom 2014)',0.090,120
    'D = 5.0 um, 30 us (Astrom 2014)',0.177,30
    'D = 5.0 um, 60 us (Astrom 2014)',0.110,60
    'D = 5.0 um, 90 us (Astrom 2014)',0.087,90
    'D = 5.0 um, 120 us (Astrom 2014)',0.074,120
    };

set(handles.ethreshpresets,'String',etv(:,1));
setappdata(handles.ethreshpresets,'data',etv);

% UIWAIT makes ea_vatsettings_horn wait for user response (see UIRESUME)

% --- Outputs from this function are returned to the command line.
function varargout = ea_vatsettings_dembek_OutputFcn(hObject, eventdata, handles)
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
vatsettings.dembek_ethresh=str2double(get(handles.ethreshtext,'String'));
vatsettings.dembek_ethreshpw=get(handles.ethreshtext,'Value');
vatsettings.dembek_pw = str2double(get(handles.pwtext,'String'));

ea_setprefs('vatsettings',vatsettings);

delete(handles.setfig);

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
    set(handles.ethreshtext,'String',num2str(data{get(hObject,'Value'),2}));
    set(handles.ethreshtext,'Value',data{get(hObject,'Value'),3});
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

% % --- Executes on selection change in pwpresets.
function pwpresets_Callback(hObject, eventdata, handles)
% hObject    handle to pwpresets (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns pwpresets contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pwpresets

data=getappdata(hObject,'data');
if isnan(data{get(hObject,'Value'),2})
    return
else
    set(handles.pwtext,'String',num2str(data{get(hObject,'Value'),2}));
%     set(handles.condwm,'String',num2str(data{get(hObject,'Value'),3}));
end
set(hObject,'value',1);


% --- Executes during object creation, after setting all properties.
function pwpresets_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pwpresets (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function pwtext_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pwtext (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function ethreshtext_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ethreshtext (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function savebutn_CreateFcn(hObject, eventdata, handles)
% hObject    handle to savebutn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called



function pwtext_Callback(hObject, eventdata, handles)
% hObject    handle to pwtext (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pwtext as text
%        str2double(get(hObject,'String')) returns contents of pwtext as a double
