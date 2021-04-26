function varargout = ea_vatsettings_butenko(varargin)
%EA_VATSETTINGS_BUTENKO MATLAB code file for ea_vatsettings_butenko.fig
%      EA_VATSETTINGS_BUTENKO, by itself, creates a new EA_VATSETTINGS_BUTENKO or raises the existing
%      singleton*.
%
%      H = EA_VATSETTINGS_BUTENKO returns the handle to a new EA_VATSETTINGS_BUTENKO or the handle to
%      the existing singleton*.
%
%      EA_VATSETTINGS_BUTENKO('Property','Value',...) creates a new EA_VATSETTINGS_BUTENKO using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to ea_vatsettings_butenko_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      EA_VATSETTINGS_BUTENKO('CALLBACK') and EA_VATSETTINGS_BUTENKO('CALLBACK',hObject,...) call the
%      local function named CALLBACK in EA_VATSETTINGS_BUTENKO.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ea_vatsettings_butenko

% Last Modified by GUIDE v2.5 15-Apr-2021 19:44:49

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ea_vatsettings_butenko_OpeningFcn, ...
                   'gui_OutputFcn',  @ea_vatsettings_butenko_OutputFcn, ...
                   'gui_LayoutFcn',  [], ...
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


% --- Executes just before ea_vatsettings_butenko is made visible.
function ea_vatsettings_butenko_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

% Choose default command line output for ea_vatsettings_butenko
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

set(handles.setfig,'name','OSS-DBS setting');

prefs=ea_prefs('');
set(handles.calcAxonActivation,'Value',prefs.machine.vatsettings.butenko_calcAxonActivation);

connectomes = ea_genmodlist([], [], [], 'dmri');

% Check multi-tract connectome
multiTractConnFolder = [ea_getconnectomebase, 'dMRI_MultiTract'];
multiTractConn = unique(cellfun(@fileparts, ea_regexpdir(multiTractConnFolder,'\.mat$'), 'Uni', 0));
multiTractConn = strrep(multiTractConn, [multiTractConnFolder,filesep], 'Multi-Tract: ')';

connectomes = [connectomes, multiTractConn];

if isempty(connectomes)
    fprintf('\n')
    warning('off', 'backtrace');
    warning('No connectome found!');
    warning('on', 'backtrace');
    set(handles.calcAxonActivation,'Value',0);
    set(handles.calcAxonActivation,'Enable','off');
    set(handles.connectomes,'String','No connectome found!');
    set(handles.connectomes,'Enable','off');
else
    set(handles.connectomes,'String',connectomes);
    connIdx = find(ismember(connectomes,prefs.machine.vatsettings.butenko_connectome));
    if ~isempty(connIdx)
        set(handles.connectomes,'Value',connIdx);
    else
        set(handles.connectomes,'Value',1);
    end
end

if numel(prefs.machine.vatsettings.butenko_axonLength) == 1 % Normal connectome
    set(handles.axonLength,'String',num2str(prefs.machine.vatsettings.butenko_axonLength));
    set(handles.fiberDiameter,'String',num2str(prefs.machine.vatsettings.butenko_fiberDiameter));
    handles.LenDSetting.Visible = 'off';
else % Multi-Tract connectome case
    set(handles.axonLength,'String','10');
    set(handles.fiberDiameter,'String','5.7');
    set(handles.axonLength,'Enable','off');
    set(handles.fiberDiameter,'Enable','off');
    handles.LenDSetting.Visible = 'on';
end

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

set(handles.ethresh,'String',num2str(prefs.machine.vatsettings.butenko_ethresh));

set(handles.interactive,'Value',prefs.machine.vatsettings.butenko_interactive);


% --- Outputs from this function are returned to the command line.
function varargout = ea_vatsettings_butenko_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in savebutn.
function savebutn_Callback(hObject, eventdata, handles)
% hObject    handle to savebutn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
prefs=ea_prefs('');

vatsettings = prefs.machine.vatsettings;
vatsettings.butenko_calcAxonActivation = get(handles.calcAxonActivation,'Value');
if strcmp(get(handles.connectomes,'Enable'), 'on')
    connectomes = get(handles.connectomes,'String');
    vatsettings.butenko_connectome = connectomes{get(handles.connectomes,'Value')};
end
if ~startsWith(vatsettings.butenko_connectome, 'Multi-Tract: ') % Values from the EditField
    vatsettings.butenko_axonLength = str2double(get(handles.axonLength,'String'));
    vatsettings.butenko_fiberDiameter = str2double(get(handles.fiberDiameter,'String'));
else % Multi-Tract connectome case, values from appdata
    vatsettings.butenko_axonLength = getappdata(handles.setfig, 'axonLength');
    vatsettings.butenko_fiberDiameter = getappdata(handles.setfig, 'fiberDiameter');
end
vatsettings.butenko_ethresh = str2double(get(handles.ethresh,'String'));
vatsettings.butenko_interactive = get(handles.interactive,'Value');
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
        % prompt = {'Enter D [um]:','Enter PW [us]:'};
        % dlgtitle = 'Specify fiber diamter D and pulse width PW';
        % dims = [1 60];
        % definput = {'3.5','60'};
        % values = inputdlg(prompt,dlgtitle,dims,definput);
        % load activation_model_3v.mat;
        % thresh = activation_model_3v(str2num(values{2}),str2num(values{1})); % get approximation
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



function fiberDiameter_Callback(hObject, eventdata, handles)
% hObject    handle to fiberDiameter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fiberDiameter as text
%        str2double(get(hObject,'String')) returns contents of fiberDiameter as a double


% --- Executes during object creation, after setting all properties.
function fiberDiameter_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fiberDiameter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function axonLength_Callback(hObject, eventdata, handles)
% hObject    handle to axonLength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of axonLength as text
%        str2double(get(hObject,'String')) returns contents of axonLength as a double


% --- Executes during object creation, after setting all properties.
function axonLength_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axonLength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in calcAxonActivation.
function calcAxonActivation_Callback(hObject, eventdata, handles)
% hObject    handle to calcAxonActivation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of calcAxonActivation


% --- Executes on selection change in connectomes.
function connectomes_Callback(hObject, eventdata, handles)
% hObject    handle to connectomes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns connectomes contents as cell array
%        contents{get(hObject,'Value')} returns selected item from connectomes
conn = hObject.String{hObject.Value};
if ~startsWith(conn, 'Multi-Tract: ')
    handles.axonLength.Enable = 'on';
    handles.fiberDiameter.Enable = 'on';
    handles.LenDSetting.Visible = 'off';
else
    handles.axonLength.Enable = 'off';
    handles.fiberDiameter.Enable = 'off';
    ea_axonact_connsetting(strrep(conn, 'Multi-Tract: ', ''), handles.setfig);
    handles.LenDSetting.Visible = 'on';
end

% --- Executes during object creation, after setting all properties.
function connectomes_CreateFcn(hObject, eventdata, handles)
% hObject    handle to connectomes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in interactive.
function interactive_Callback(hObject, eventdata, handles)
% hObject    handle to interactive (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of interactive


% --- Executes on button press in LenDSetting.
function LenDSetting_Callback(hObject, eventdata, handles)
% hObject    handle to LenDSetting (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
conn = handles.connectomes.String{handles.connectomes.Value};
ea_axonact_connsetting(strrep(conn, 'Multi-Tract: ', ''), handles.setfig);
