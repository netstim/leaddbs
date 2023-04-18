function varargout = lead_mapper(varargin)
% LEAD_MAPPER MATLAB code for lead_mapper.fig
%      LEAD_MAPPER, by itself, creates a new LEAD_MAPPER or raises the existing
%      singleton*.
%
%      H = LEAD_MAPPER returns the handle to a new LEAD_MAPPER or the handle to
%      the existing singleton*.
%
%      LEAD_MAPPER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in LEAD_MAPPER.M with the given input arguments.
%
%      LEAD_MAPPER('Property','Value',...) creates a new LEAD_MAPPER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before lead_mapper_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to lead_mapper_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help lead_mapper

% Last Modified by GUIDE v2.5 02-Mar-2023 19:11:00

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @lead_mapper_OpeningFcn, ...
                   'gui_OutputFcn',  @lead_mapper_OutputFcn, ...
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


% --- Executes just before lead_mapper is made visible.
function lead_mapper_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to lead_mapper (see VARARGIN)

handles.prod = 'mapper';
handles.callingfunction = 'lead_mapper';

earoot=ea_getearoot;
im=imread([earoot,'icons',filesep,'logo_lead_connectome_mapper.png']);
image(im);
axes(handles.logoaxes);
axis off;
axis equal;
set(handles.leadfigure,'name','Lead Connectome Mapper','color','w');
% homedir=ea_gethome;
%setappdata(handles.leadfigure,'uipatdir',{homedir(1:end-1)});

% add recentpatients patients...
ea_initrecent(handles, 'patients');

ea_processguiargs(handles,varargin)

ea_menu_initmenu(handles,{'cluster','prefs','transfer','vats'},ea_prefs);


[mdl,sf]=ea_genmodlist;
ea_updatemodpopups(mdl,sf,handles)

set(handles.versiontxt,'String',['v',ea_getvsn('local')]);

ea_bind_dragndrop(handles.leadfigure, ...
    @(obj,evt) DropFcn(obj,evt,handles), ...
    @(obj,evt) DropFcn(obj,evt,handles));

ea_ListBoxRenderer(handles.fiberspopup);
ea_ListBoxRenderer(handles.fmripopup);

% Choose default command line output for lead_mapper
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% Disable buttons for standalone app
if isdeployed
    set(handles.exportcode,'Enable','off');
end

% UIWAIT makes lead_mapper wait for user response (see UIRESUME)
% uiwait(handles.leadfigure);


% --- Drag and drop callback to load patdir.
function DropFcn(~, event, handles)

switch event.DropType
    case 'file'
        patdir = event.Data;
    case 'string'
        patdir = {event.Data};
end

nonexist = cellfun(@(x) ~exist(x, 'dir'), patdir);
if any(nonexist)
    fprintf('\nExcluded non-existent/invalid folder:\n');
    cellfun(@disp, patdir(nonexist));
    fprintf('\n');
    patdir(nonexist) = [];
end

ea_busyaction('on',handles.leadfigure,'mapper');
if ~isempty(patdir)
    ea_load_pts(handles, patdir);
end
ea_busyaction('off',handles.leadfigure,'mapper');


% --- Outputs from this function are returned to the command line.
function varargout = lead_mapper_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in fiberspopup.
function fiberspopup_Callback(hObject, eventdata, handles)
% hObject    handle to fiberspopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns fiberspopup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from fiberspopup


% --- Executes during object creation, after setting all properties.
function fiberspopup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fiberspopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu5.
function popupmenu5_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu5 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu5


% --- Executes during object creation, after setting all properties.
function popupmenu5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox5.
function checkbox5_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox5


% --- Executes on button press in dostructural.
function dostructural_Callback(hObject, eventdata, handles)
% hObject    handle to dostructural (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of dostructural


% --- Executes on button press in seedbutton.
function seedbutton_Callback(hObject, eventdata, handles)
% hObject    handle to seedbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

seeddef = get(handles.seeddefpopup,'String');
switch seeddef{get(handles.seeddefpopup,'Value')}
    case 'Manually choose seeds'
        [seeds,path] = uigetfile({'*'},'Please choose seed definition(s)...','MultiSelect','on');
    case 'Manually choose parcellation'
        [seeds,path] = uigetfile({'*'},'Please choose parcellation...',ea_space([],'labeling'),'MultiSelect','off');
end

if ischar(path) % path is 0 if the user clicks Cancel or close the window
    if iscell(seeds)
        set(hObject,'String',['Multiple (',num2str(length(seeds)),')']);
    elseif ischar(seeds)
        set(hObject,'String',seeds);
        seeds = {seeds};
    end

    for s = 1:length(seeds)
        seeds{s} = fullfile(path,seeds{s});
    end

    setappdata(hObject,'seeds',seeds);
end


% --- Executes on button press in run_button.
function run_button_Callback(hObject, eventdata, handles)
% hObject    handle to run_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

leadfigure = handles.leadfigure;
ea_busyaction('on',leadfigure,'mapper');

options = ea_handles2options(handles);
options.uivatdirs = getappdata(handles.leadfigure,'uipatdir');
options.uipatdirs = {''};

options.leadprod = 'mapper';

ea_run('run',options);

ea_busyaction('off',leadfigure,'mapper');


% --- Executes on button press in exportcode.
function exportcode_Callback(hObject, eventdata, handles)
% hObject    handle to exportcode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

leadfigure = handles.leadfigure;
ea_busyaction('on',leadfigure,'mapper');

options = ea_handles2options(handles);
options.uivatdirs = getappdata(handles.leadfigure,'uipatdir');
options.uipatdirs = {''};

options.leadprod = 'mapper';

ea_run('export', options);

ea_busyaction('off',leadfigure,'mapper');


% --- Executes on selection change in command.
function command_Callback(hObject, eventdata, handles)
% hObject    handle to command (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns command contents as cell array
%        contents{get(hObject,'Value')} returns selected item from command
conn = handles.fmripopup.String{handles.fmripopup.Value};
cmd = handles.command.String{handles.command.Value};
if contains(cmd, 'matrix', 'IgnoreCase', true) && ~contains(conn, 'matrix', 'IgnoreCase', true)
    handles.exportgmtc.Visible = 'on';
else
    handles.exportgmtc.Visible = 'off';
end


% --- Executes during object creation, after setting all properties.
function command_CreateFcn(hObject, eventdata, handles)
% hObject    handle to command (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in dofunctional.
function dofunctional_Callback(hObject, eventdata, handles)
% hObject    handle to dofunctional (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of dofunctional
conn = handles.fmripopup.String{handles.fmripopup.Value};
cmd = handles.command.String{handles.command.Value};
if contains(cmd, 'matrix', 'IgnoreCase', true) && ~contains(conn, 'matrix', 'IgnoreCase', true)
    handles.exportgmtc.Visible = 'on';
else
    handles.exportgmtc.Visible = 'off';
end


% --- Executes on selection change in fmripopup.
function fmripopup_Callback(hObject, eventdata, handles)
% hObject    handle to fmripopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns fmripopup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from fmripopup
conn = handles.fmripopup.String{handles.fmripopup.Value};
cmd = handles.command.String{handles.command.Value};
if contains(cmd, 'matrix', 'IgnoreCase', true) && ~contains(conn, 'matrix', 'IgnoreCase', true)
    handles.exportgmtc.Visible = 'on';
else
    handles.exportgmtc.Visible = 'off';
end


% --- Executes during object creation, after setting all properties.
function fmripopup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fmripopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in odirbutton.
function odirbutton_Callback(hObject, eventdata, handles)
% hObject    handle to odirbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

seeds=getappdata(handles.seedbutton,'seeds');
if ~isempty(seeds) % seeds defined already
    seedbase=fileparts(seeds{1});
else
    seedbase='';
end
odir=uigetdir(seedbase,'Choose output location');
if ischar(odir)
    setappdata(hObject,'odir',[odir,filesep]);
    set(hObject,'String',odir);
end


% --- Executes on selection change in strucexportspace.
function strucexportspace_Callback(hObject, eventdata, handles)
% hObject    handle to strucexportspace (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns strucexportspace contents as cell array
%        contents{get(hObject,'Value')} returns selected item from strucexportspace


% --- Executes during object creation, after setting all properties.
function strucexportspace_CreateFcn(hObject, eventdata, handles)
% hObject    handle to strucexportspace (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in omaskbutton.
function omaskbutton_Callback(hObject, eventdata, handles)
% hObject    handle to omaskbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[omask, path] = uigetfile({'*.nii';'*.nii.gz'},'Choose output location');
if ischar(path)
    setappdata(hObject,'omask',fullfile(path,omask));
    set(hObject,'String',omask);
end


% --- Executes on button press in patdir_choosebox.
function patdir_choosebox_Callback(hObject, eventdata, handles)
% hObject    handle to patdir_choosebox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ea_busyaction('on',handles.leadfigure,'mapper');
options.prefs = ea_prefs;
ea_getpatients(options,handles);
ea_busyaction('off',handles.leadfigure,'mapper');

% --- Executes on selection change in recentpatients.
function recentpatients_Callback(hObject, eventdata, handles)
% hObject    handle to recentpatients (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns recentpatients contents as cell array
%        contents{get(hObject,'Value')} returns selected item from recentpatients
ea_busyaction('on',handles.leadfigure,'mapper');
ea_recentcallback(handles, 'patients');
ea_busyaction('off',handles.leadfigure,'mapper');


% --- Executes during object creation, after setting all properties.
function recentpatients_CreateFcn(hObject, eventdata, handles)
% hObject    handle to recentpatients (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in seeddefpopup.
function seeddefpopup_Callback(hObject, eventdata, handles)
% hObject    handle to seeddefpopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns seeddefpopup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from seeddefpopup

str=get(hObject,'String');
if iscell(str)
    str=str{get(hObject,'Value')};
end
if strcmp(str, 'Manually choose seeds')
   set(handles.seedbutton,'enable','on');
   set(handles.seedbutton,'String','Choose seeds...');
elseif strcmp(str, 'Manually choose parcellation')
   set(handles.seedbutton,'enable','on');
   set(handles.seedbutton,'String','Choose parcellation...');
else
   set(handles.seedbutton,'enable','off');
   set(handles.seedbutton,'String','Choose seeds/parcellation...');
end

% --- Executes during object creation, after setting all properties.
function seeddefpopup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to seeddefpopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in openpatientdir.
function openpatientdir_Callback(hObject, eventdata, handles)
% hObject    handle to openpatientdir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ea_openpatdir(handles);


% --- Executes on button press in exportgmtc.
function exportgmtc_Callback(hObject, eventdata, handles)
% hObject    handle to exportgmtc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of exportgmtc
