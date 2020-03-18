function varargout = ea_vatsettings_fastfield(varargin)
% EA_VATSETTINGS_FASTFIELD MATLAB code for ea_vatsettings_fastfield.fig
%      EA_VATSETTINGS_FASTFIELD, by itself, creates a new EA_VATSETTINGS_FASTFIELD or raises the existing
%      singleton*.
%
%      H = EA_VATSETTINGS_FASTFIELD returns the handle to a new EA_VATSETTINGS_FASTFIELD or the handle to
%      the existing singleton*.
%
%      EA_VATSETTINGS_FASTFIELD('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in EA_VATSETTINGS_FASTFIELD.M with the given input arguments.
%
%      EA_VATSETTINGS_FASTFIELD('Property','Value',...) creates a new EA_VATSETTINGS_FASTFIELD or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ea_vatsettings_fastfield_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ea_vatsettings_fastfield_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ea_vatsettings_fastfield

% Last Modified by GUIDE v2.5 11-Mar-2020 17:00:53

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ea_vatsettings_fastfield_OpeningFcn, ...
                   'gui_OutputFcn',  @ea_vatsettings_fastfield_OutputFcn, ...
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


% --- Executes just before ea_vatsettings_fastfield is made visible.
function ea_vatsettings_fastfield_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ea_vatsettings_fastfield (see VARARGIN)

% Choose default command line output for ea_vatsettings_fastfield
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);


set(handles.setfig,'name','FEM-based VAT model setting');

prefs=ea_prefs('');

set(handles.condb,'String',num2str(prefs.machine.vatsettings.fastfield_cb));
set(handles.ethresh,'String',num2str(prefs.machine.vatsettings.fastfield_ethresh));  


 ea_fillpresetpopups(handles);


% UIWAIT makes ea_vatsettings_fastfield wait for user response (see UIRESUME)
% uiwait(handles.setfig);


function ea_fillpresetpopups(handles)


etv={'E-Field Threshold Presets:',nan
    'Approximation by D [um] and PW [us] (Proverbio & Husch 2019)',nan
    'General Heuristic (e.g. Hemm 2005. Vasques 2009, Astrom 2009, Horn 2017)',0.2
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


% --- Outputs from this function are returned to the command line.
function varargout = ea_vatsettings_fastfield_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function condb_Callback(hObject, eventdata, handles)
% hObject    handle to condb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of condb as text
%        str2double(get(hObject,'String')) returns contents of condb as a double


% --- Executes during object creation, after setting all properties.
function condb_CreateFcn(hObject, eventdata, handles)
% hObject    handle to condb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



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
thresh=data{get(hObject,'Value'),2};
if isnan(thresh)
    if(strcmp(data{get(hObject,'Value'),1}, ...
            'Approximation by D [um] and PW [us] (Proverbio & Husch 2019)'))
        f = approxonGui;
        uiwait(f); % setting thresh
        thresh = getappdata(f, 'thresh');
        close(f);
    else
        return
    end
end
set(handles.ethresh,'String',num2str(thresh));

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


% --- Executes on button press in savebutn.
function savebutn_Callback(hObject, eventdata, handles, varargin)
% hObject    handle to savebutn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
prefs=ea_prefs('');

vatsettings=prefs.machine.vatsettings;
vatsettings.fastfield_cb=str2double(get(handles.condb,'String'));
vatsettings.fastfield_ethresh=str2double(get(handles.ethresh,'String'));
ea_setprefs('vatsettings',vatsettings);

delete(handles.setfig);
