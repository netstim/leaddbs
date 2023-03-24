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

% Last Modified by GUIDE v2.5 22-Jan-2020 08:56:43

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
set(handles.useatlas,'Value',prefs.machine.vatsettings.horn_useatlas);
ea_refreshgmwm(handles);
options=ea_defaultoptions;
options.prefs.machine.defaultatlas=prefs.machine.vatsettings.horn_atlasset;
ea_listatlassets(options,handles,1);

set(handles.removeElectrode,'Value',prefs.machine.vatsettings.horn_removeElectrode);

ea_fillpresetpopups(handles);


function ea_fillpresetpopups(handles)

condv={'Conductivity Presets:',nan,nan
    'General Heuristic (e.g. Buzsaki 2006)',0.33,0.14
    'Frequency Adapted (10 kHz adj. cf. itis.ethz.ch)',0.1365,0.09104};  % adjusted for 1 kOhm at 10 kHz for 0.1 mm GM encap (no EDL)
set(handles.condpresets,'String',condv(:,1));
setappdata(handles.condpresets,'data',condv);


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
    'D = 4.0 um, 60 us (Astrom 2014)',0.150
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
vatsettings.horn_useatlas=get(handles.useatlas,'Value');
vatsettings.horn_atlasset=get(handles.atlassetpopup,'String');
vatsettings.horn_atlasset=vatsettings.horn_atlasset{get(handles.atlassetpopup,'Value')};
vatsettings.horn_removeElectrode=get(handles.removeElectrode,'Value');
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
thresh=data{get(hObject,'Value'),2};
if isnan(thresh)
    if(strcmp(data{get(hObject,'Value'),1}, ...
            'Approximation by D [um] and PW [us] (Proverbio & Husch 2019)'))
        % Prompt an input dialog
%         prompt = {'Enter D [um]:','Enter PW [us]:'};
%         dlgtitle = 'Specify fiber diamter D and pulse width PW';
%         dims = [1 60];
%         definput = {'3.5','60'};
%         values = inputdlg(prompt,dlgtitle,dims,definput);
%         load activation_model_3v.mat;
%         thresh = activation_model_3v(str2num(values{2}),str2num(values{1})); % get approximation
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



function condgm_Callback(hObject, eventdata, handles)
% hObject    handle to condgm_txt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of condgm_txt as text
%        str2double(get(hObject,'String')) returns contents of condgm_txt as a double


% --- Executes during object creation, after setting all properties.
function condgm_CreateFcn(hObject, eventdata, handles)
% hObject    handle to condgm_txt (see GCBO)
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
% hObject    handle to condwm_txt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of condwm_txt as text
%        str2double(get(hObject,'String')) returns contents of condwm_txt as a double


% --- Executes during object creation, after setting all properties.
function condwm_CreateFcn(hObject, eventdata, handles)
% hObject    handle to condwm_txt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in useatlas.
function useatlas_Callback(hObject, eventdata, handles)
% hObject    handle to useatlas (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of useatlas
ea_refreshgmwm(handles);

function ea_refreshgmwm(handles)
switch get(handles.useatlas,'Value')
    case 1
        set(handles.condgm_txt,'visible','on');
        set(handles.condgm,'visible','on');
        set(handles.smgm_txt,'visible','on');
        set(handles.atlassetpopup,'visible','on');
        set(handles.condwm_txt,'string','Conductivity, WM:');
    case 0
        set(handles.condgm_txt,'visible','off');
        set(handles.condgm,'visible','off');
        set(handles.smgm_txt,'visible','off');
        set(handles.atlassetpopup,'visible','off');
        set(handles.condwm_txt,'string','Conductivity:');
end


% --- Executes on selection change in atlassetpopup.
function atlassetpopup_Callback(hObject, eventdata, handles)
% hObject    handle to atlassetpopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns atlassetpopup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from atlassetpopup


% --- Executes during object creation, after setting all properties.
function atlassetpopup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to atlassetpopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in removeElectrode.
function removeElectrode_Callback(hObject, eventdata, handles)
% hObject    handle to removeElectrode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of removeElectrode
ea_setprefs('vatsettings.horn_removeElectrode',get(handles.removeElectrode,'Value'));
