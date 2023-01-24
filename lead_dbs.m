function varargout = lead_dbs(varargin)
% LEAD_DBS M-file for lead_dbs.fig
%      LEAD_DBS, by itself, creates a new LEAD_DBS or raises the existing
%      singleton*.
%
%      H = LEAD_DBS returns the handle to a new LEAD_DBS or the handle to
%      the existing singleton*.
%
%      LEAD_DBS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in LEAD_DBS.M with the given input arguments.
%
%      LEAD_DBS('Property','Value',...) creates a new LEAD_DBS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before lead_dbs_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to lead_dbs_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help lead_dbs

% Last Modified by GUIDE v2.5 16-Sep-2022 19:45:41

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @lead_dbs_OpeningFcn, ...
    'gui_OutputFcn',  @lead_dbs_OutputFcn, ...
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


% --- Executes just before lead_dbs is made visible.
function lead_dbs_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to lead_dbs (see VARARGIN)

earoot=ea_getearoot;

handles.prod = 'dbs';

% add recent patients...
ea_initrecent(handles, 'patients');

set(handles.vizspacepopup,'String',{[ea_underscore2space(ea_getspace), ' Space'];'Native Patient Space'});

ea_dispbn;

set(handles.leadfigure,'name','Welcome to LEAD-DBS');

spacedef=ea_getspacedef;
if isfield(spacedef,'guidef')
    set(handles.targetpopup,'String',[spacedef.guidef.entrypoints,{'Manual'}, {'Auto'}]);
end

options.prefs=ea_prefs('');

ea_init_coregmrpopup(handles, options.prefs.mrcoreg.default);
ea_init_coregctpopup(handles, options.prefs.ctcoreg.default);

% load atlassets
ea_listatlassets(options,handles,1);

set(handles.normalize_checkbox,'Value',0);

set(hObject,'Color',[1 1 1]);
set(handles.versiontxt,'String',['v',ea_getvsn('local')]);

set(0,'CurrentFigure',handles.leadfigure);
im=imread([earoot,'icons',filesep,'logo_lead_dbs.png']);

image(im);
axis off;
axis equal;

% get electrode model specs and place in popup
set(handles.electrode_model_popup,'String',ea_resolve_elspec);

options.earoot=earoot;

% Initialize norm methods popupmenu
ea_init_normpopup(handles, options.prefs.normalize.default);

ea_processguiargs(handles,varargin)

%% add tools menu
ea_menu_initmenu(handles,{'import','acpc','export','applynorm','leador','dbs','cluster','prefs','vatcon','transfer','checkregfigs','space','surfice','methods'},options.prefs);

ea_firstrun(handles,options);
ea_getui(handles);

ea_bind_dragndrop(handles.leadfigure, ...
    @(obj,evt) DropFcn(obj,evt,handles), ...
    @(obj,evt) DropFcn(obj,evt,handles));

ea_ListBoxRenderer(handles.electrode_model_popup);
ea_ListBoxRenderer(handles.normmethod);
ea_ListBoxRenderer(handles.scrfmask);
ea_ListBoxRenderer(handles.atlassetpopup);

% Choose default command line output for lead_dbs
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% Disable buttons for standalone app
if isdeployed
    set(handles.updatebutn,'Enable','off');
end

% UIWAIT makes lead_dbs wait for user response (see UIRESUME)
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

ea_busyaction('on',handles.leadfigure,'dbs');
if ~isempty(patdir)
    ea_load_pts(handles, patdir);

    if isfield(handles,'atlassetpopup')
        atlasset=get(handles.atlassetpopup,'String');
        atlasset=atlasset{get(handles.atlassetpopup,'Value')};
        options.prefs=ea_prefs('');
        ea_listatlassets(options,handles,get(handles.vizspacepopup,'Value'),atlasset);
    end
end
ea_busyaction('off',handles.leadfigure,'dbs');


% --- Outputs from this function are returned to the command line.
function varargout = lead_dbs_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
try
    varargout{1} = handles.output;
end

% --- Executes on button press in run_button.
function run_button_Callback(hObject, eventdata, handles)
% hObject    handle to run_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

ea_busyaction('on', handles.leadfigure, 'dbs');

options = ea_handles2options(handles);
options.uipatdirs = getappdata(handles.leadfigure,'uipatdir');

options.leadprod = 'dbs';

setappdata(handles.leadfigure,'handles',handles);
options.leadfigure=handles.leadfigure;

ea_run('run',options);

ea_busyaction('off', handles.leadfigure, 'dbs');


function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double
ea_storeui(handles);


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in electrode_model_popup.
function electrode_model_popup_Callback(hObject, eventdata, handles)
% hObject    handle to electrode_model_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns electrode_model_popup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from electrode_model_popup
ea_storeui(handles);


% --- Executes during object creation, after setting all properties.
function electrode_model_popup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to electrode_model_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function endtol_txt_Callback(hObject, eventdata, handles)
% hObject    handle to endtol_txt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of endtol_txt as text
%        str2double(get(hObject,'String')) returns contents of endtol_txt as a double
ea_storeui(handles);


% --- Executes during object creation, after setting all properties.
function endtol_txt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to endtol_txt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function distance_txt_Callback(hObject, eventdata, handles)
% hObject    handle to distance_txt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of distance_txt as text
%        str2double(get(hObject,'String')) returns contents of distance_txt as a double
ea_storeui(handles);


% --- Executes during object creation, after setting all properties.
function distance_txt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to distance_txt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function axis_std_txt_Callback(hObject, eventdata, handles)
% hObject    handle to axis_std_txt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of axis_std_txt as text
%        str2double(get(hObject,'String')) returns contents of axis_std_txt as a double
ea_storeui(handles);


% --- Executes during object creation, after setting all properties.
function axis_std_txt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axis_std_txt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function cor_std_txt_Callback(hObject, eventdata, handles)
% hObject    handle to cor_std_txt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of cor_std_txt as text
%        str2double(get(hObject,'String')) returns contents of cor_std_txt as a double
ea_storeui(handles);


% --- Executes during object creation, after setting all properties.
function cor_std_txt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cor_std_txt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function axiscontrast_txt_Callback(hObject, eventdata, handles)
% hObject    handle to axiscontrast_txt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of axiscontrast_txt as text
%        str2double(get(hObject,'String')) returns contents of axiscontrast_txt as a double
ea_storeui(handles);


% --- Executes during object creation, after setting all properties.
function axiscontrast_txt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axiscontrast_txt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in slowdemo_checkbox.
function slowdemo_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to slowdemo_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of slowdemo_checkbox
ea_storeui(handles);



function z_contrast_txt_Callback(hObject, eventdata, handles)
% hObject    handle to z_contrast_txt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of z_contrast_txt as text
%        str2double(get(hObject,'String')) returns contents of z_contrast_txt as a double
ea_storeui(handles);


% --- Executes during object creation, after setting all properties.
function z_contrast_txt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to z_contrast_txt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in prolong_electrode_checkbox.
function prolong_electrode_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to prolong_electrode_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of prolong_electrode_checkbox
ea_storeui(handles);


% --- Executes on button press in render_checkbox.
function render_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to render_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of render_checkbox
ea_storeui(handles);


% --- Executes on button press in showatlases_checkbox.
function showatlases_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to showatlases_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of showatlases_checkbox
ea_storeui(handles);


% --- Executes on button press in showfibers_checkbox.
function showfibers_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to showfibers_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of showfibers_checkbox
ea_storeui(handles);


% --- Executes on button press in showconnectivities_checkbox.
function showconnectivities_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to showconnectivities_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of showconnectivities_checkbox
ea_storeui(handles);


% --- Executes on button press in writeout2d_checkbox.
function writeout2d_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to writeout2d_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of writeout2d_checkbox
ea_storeui(handles);


% --- Executes on button press in showatlases2d_checkbox.
function showatlases2d_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to showatlases2d_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of showatlases2d_checkbox
ea_storeui(handles);


% --- Executes on button press in manualheight_checkbox.
function manualheight_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to manualheight_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of manualheight_checkbox
ea_storeui(handles);


% --- Executes on button press in patdir_choosebox.
function patdir_choosebox_Callback(hObject, eventdata, handles)
% hObject    handle to patdir_choosebox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ea_busyaction('on',handles.leadfigure,'dbs');
options.prefs = ea_prefs;
ea_getpatients(options,handles);
ea_busyaction('off',handles.leadfigure,'dbs');


% --- Executes on button press in left_checkbox.
function left_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to left_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of left_checkbox
ea_storeui(handles);


% --- Executes on button press in right_checkbox.
function right_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to right_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of right_checkbox
ea_storeui(handles);


% --- Executes on button press in normalize_checkbox.
function normalize_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to normalize_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of normalize_checkbox
ea_storeui(handles);



function refinesteps_txt_Callback(~, eventdata, handles)
% hObject    handle to refinesteps_txt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of refinesteps_txt as text
%        str2double(get(hObject,'String')) returns contents of refinesteps_txt as a double
ea_storeui(handles);


% --- Executes during object creation, after setting all properties.
function refinesteps_txt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to refinesteps_txt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end








function maskwindow_txt_Callback(hObject, eventdata, handles)
% hObject    handle to maskwindow_txt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of maskwindow_txt as text
%        str2double(get(hObject,'String')) returns contents of maskwindow_txt as a double
ea_storeui(handles);


% --- Executes during object creation, after setting all properties.
function maskwindow_txt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to maskwindow_txt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in doreconstruction_checkbox.
function doreconstruction_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to doreconstruction_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of doreconstruction_checkbox

if get(handles.reconmethod,'Value')==1 && get(hObject,'Value')==1
    set(handles.maskwindow_txt,'Enable','on');
    set(handles.targetpopup,'Enable','on');
else
    set(handles.maskwindow_txt,'Enable','off');
    set(handles.targetpopup,'Enable','off');
end

ea_storeui(handles);


% --- Executes on selection change in MRCT.
function MRCT_Callback(hObject, eventdata, handles)
% hObject    handle to MRCT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns MRCT contents as cell array
%        contents{get(hObject,'Value')} returns selected item from MRCT

ea_switchctmr(handles, get(hObject,'Value'));

ea_storeui(handles);


% --- Executes during object creation, after setting all properties.
function MRCT_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MRCT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in stimcheckbox.
function stimcheckbox_Callback(hObject, eventdata, handles)
% hObject    handle to stimcheckbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of stimcheckbox
ea_storeui(handles);


% --- Executes on button press in cg25check.
function cg25check_Callback(hObject, eventdata, handles)
% hObject    handle to cg25check (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cg25check
ea_storeui(handles);


% --- Executes on selection change in normmethod.
function normmethod_Callback(hObject, eventdata, handles)
% hObject    handle to normmethod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns normmethod contents as cell array
%        contents{get(hObject,'Value')} returns selected item from normmethod
ea_storeui(handles);
ea_normsettings(handles);


% --- Executes during object creation, after setting all properties.
function normmethod_CreateFcn(hObject, eventdata, handles)
% hObject    handle to normmethod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in axispopup.
function axispopup_Callback(hObject, eventdata, handles)
% hObject    handle to axispopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns axispopup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from axispopup
ea_storeui(handles);


% --- Executes during object creation, after setting all properties.
function axispopup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axispopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on selection change in zpopup.
function zpopup_Callback(hObject, eventdata, handles)
% hObject    handle to zpopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns zpopup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from zpopup
ea_storeui(handles);


% --- Executes during object creation, after setting all properties.
function zpopup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to zpopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on selection change in targetpopup.
function targetpopup_Callback(hObject, eventdata, handles)
% hObject    handle to targetpopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns targetpopup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from targetpopup
ea_storeui(handles);


% --- Executes during object creation, after setting all properties.
function targetpopup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to targetpopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in eepush.
function eepush_Callback(hObject, eventdata, handles)
% hObject    handle to eepush (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ea_dispbn('ee');



% --- Executes on selection change in atlassetpopup.
function atlassetpopup_Callback(hObject, eventdata, handles)
% hObject    handle to atlassetpopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns atlassetpopup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from atlassetpopup
ea_storeui(handles);


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


% --- Executes on button press in normcheck.
function normcheck_Callback(hObject, eventdata, handles)
% hObject    handle to normcheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of normcheck
ea_storeui(handles);


% --- Executes on button press in cmappushbutton.
function cmappushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to cmappushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
colormapeditor;
ea_storeui(handles);



% --- Executes on button press in canatlcheck.
function canatlcheck_Callback(hObject, eventdata, handles)
% hObject    handle to canatlcheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of canatlcheck
ea_storeui(handles);


% --- Executes on button press in patatlcheck.
function patatlcheck_Callback(hObject, eventdata, handles)
% hObject    handle to patatlcheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of patatlcheck
ea_storeui(handles);


% --- Executes on button press in normptatlascheck.
function normptatlascheck_Callback(hObject, eventdata, handles)
% hObject    handle to normptatlascheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of normptatlascheck
ea_storeui(handles);


% --- Executes on button press in updatebutn.
function updatebutn_Callback(hObject, eventdata, handles)
% hObject    handle to updatebutn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ea_update;


% --- Executes on button press in exportservercheck.
function exportservercheck_Callback(hObject, eventdata, handles)
% hObject    handle to exportservercheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of exportservercheck
ea_storeui(handles);


% --- Executes on button press in coregct_checkbox.
function coregct_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to coregct_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of coregct_checkbox
ea_storeui(handles);


% --- Executes on selection change in coregctmethod.
function coregctmethod_Callback(hObject, eventdata, handles)
% hObject    handle to coregctmethod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns coregctmethod contents as cell array
%        contents{get(hObject,'Value')} returns selected item from coregctmethod
methods=getappdata(handles.leadfigure,'coregctmethod');
ea_storeui(handles);


% --- Executes during object creation, after setting all properties.
function coregctmethod_CreateFcn(hObject, eventdata, handles)
% hObject    handle to coregctmethod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
ea_storeui(handles);


% --- Executes on button press in coregctcheck.
function coregctcheck_Callback(hObject, eventdata, handles)
% hObject    handle to coregctcheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of coregctcheck
ea_storeui(handles);


% --- Executes on button press in ft_checkbox.
function ft_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to ft_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ft_checkbox
ea_storeui(handles);

% --- Executes on selection change in ftmethod.
function ftmethod_Callback(hObject, eventdata, handles)
% hObject    handle to ftmethod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns ftmethod contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ftmethod
ea_storeui(handles);

% --- Executes during object creation, after setting all properties.
function ftmethod_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ftmethod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in parcellation_atlas.
function parcellation_atlas_Callback(hObject, eventdata, handles)
% hObject    handle to parcellation_atlas (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns parcellation_atlas contents as cell array
%        contents{get(hObject,'Value')} returns selected item from parcellation_atlas
ea_storeui(handles);


% --- Executes during object creation, after setting all properties.
function parcellation_atlas_CreateFcn(hObject, eventdata, handles)
% hObject    handle to parcellation_atlas (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in include_lead_connectome_subroutine.
function include_lead_connectome_subroutine_Callback(hObject, eventdata, handles)
% hObject    handle to include_lead_connectome_subroutine (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of include_lead_connectome_subroutine


% --- Executes on button press in openleadconnectome.
function openleadconnectome_Callback(hObject, eventdata, handles)
% hObject    handle to openleadconnectome (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
lead_connectome('dependent');


% --- Executes on button press in specify2dwrite.
function specify2dwrite_Callback(hObject, eventdata, handles)
% hObject    handle to specify2dwrite (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
options = ea_handles2options(handles);
if ~isempty(getappdata(handles.leadfigure,'uipatdir'))
    bids = getappdata(handles.leadfigure, 'bids');
    subjId = getappdata(handles.leadfigure, 'subjId');
    options.subj = bids.getSubj(subjId{1}, options.modality);
end
ea_spec2dwrite(options);


% --- Executes on button press in importfs.
function importfs_Callback(hObject, eventdata, handles)
% hObject    handle to importfs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ea_importfs(handles)


% --- Executes on button press in viewmanual.
function viewmanual_Callback(hObject, eventdata, handles)
% hObject    handle to viewmanual (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
web('https://www.lead-dbs.org/?page_id=71', '-browser')


% --- Executes on selection change in vizspacepopup.
function vizspacepopup_Callback(hObject, eventdata, handles)
% hObject    handle to vizspacepopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns vizspacepopup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from vizspacepopup

% if get(hObject,'Value')==2
%    set(handles.writeout2d_checkbox,'Enable','off');
%    set(handles.writeout2d_checkbox,'Value',0);
% else
%    set(handles.writeout2d_checkbox,'Enable','on');
%    %set(handles.writeout2d_checkbox,'Value',1);
% end
atlasset=get(handles.atlassetpopup,'String');
if get(handles.atlassetpopup,'Value')>length(atlasset)
    set(handles.atlassetpopup,'Value',length(atlasset));
end
atlasset=atlasset{get(handles.atlassetpopup,'Value')};
options.prefs=ea_prefs('');
ea_listatlassets(options,handles,get(handles.vizspacepopup,'Value'),atlasset);


% --- Executes during object creation, after setting all properties.
function vizspacepopup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to vizspacepopup (see GCBO)
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



% --- Executes on selection change in coregmrmethod.
function coregmrmethod_Callback(hObject, eventdata, handles)
% hObject    handle to coregmrmethod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns coregmrmethod contents as cell array
%        contents{get(hObject,'Value')} returns selected item from coregmrmethod


% --- Executes during object creation, after setting all properties.
function coregmrmethod_CreateFcn(hObject, eventdata, handles)
% hObject    handle to coregmrmethod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in exportcode.
function exportcode_Callback(hObject, eventdata, handles)
% hObject    handle to exportcode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

ea_busyaction('on', handles.leadfigure, 'dbs');

options = ea_handles2options(handles);
options.uipatdirs = getappdata(handles.leadfigure,'uipatdir');

options.leadprod = 'dbs';

ea_run('export',options);

ea_busyaction('off', handles.leadfigure, 'dbs');


% --- Executes during object creation, after setting all properties.
function leadfigure_CreateFcn(hObject, eventdata, handles)
% hObject    handle to leadfigure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
label='lead-dbs.org';
url='https://www.lead-dbs.org/';
position=[65,hObject.Position(4)-60,85,16];
ea_hyperlink_label(label, url, position);


% --- Executes on selection change in recent.
function recent_Callback(hObject, eventdata, handles)
% hObject    handle to recent (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns recent contents as cell array
%        contents{get(hObject,'Value')} returns selected item from recent
ea_busyaction('on',handles.leadfigure,'dbs');
ea_recentcallback(handles, 'patients');
ea_busyaction('off',handles.leadfigure,'dbs');

% --- Executes during object creation, after setting all properties.
function recent_CreateFcn(hObject, eventdata, handles)
% hObject    handle to recent (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function normalize_checkbox_ButtonDownFcn(hObject,eventdata,handles)
ea_gethelp(get(handles.leadfigure,'SelectionType'),hObject);


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over normcheck.
function normcheck_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to normcheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ea_gethelp(get(handles.leadfigure,'SelectionType'),hObject);


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over doreconstruction_checkbox.
function doreconstruction_checkbox_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to doreconstruction_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ea_gethelp(get(handles.leadfigure,'SelectionType'),hObject);


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over manualheight_checkbox.
function manualheight_checkbox_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to manualheight_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ea_gethelp(get(handles.leadfigure,'SelectionType'),hObject);


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over include_lead_connectome_subroutine.
function include_lead_connectome_subroutine_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to include_lead_connectome_subroutine (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ea_gethelp(get(handles.leadfigure,'SelectionType'),hObject);


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over coregct_checkbox.
function coregct_checkbox_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to coregct_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ea_gethelp(get(handles.leadfigure,'SelectionType'),hObject);


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over coregctcheck.
function coregctcheck_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to coregctcheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ea_gethelp(get(handles.leadfigure,'SelectionType'),hObject);


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over writeout2d_checkbox.
function writeout2d_checkbox_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to writeout2d_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ea_gethelp(get(handles.leadfigure,'SelectionType'),hObject);


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over render_checkbox.
function render_checkbox_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to render_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ea_gethelp(get(handles.leadfigure,'SelectionType'),hObject);


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over exportservercheck.
function exportservercheck_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to exportservercheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ea_gethelp(get(handles.leadfigure,'SelectionType'),hObject);


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over run_button.
function run_button_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to run_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ea_gethelp(get(handles.leadfigure,'SelectionType'),hObject);


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over exportcode.
function exportcode_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to exportcode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ea_gethelp(get(handles.leadfigure,'SelectionType'),hObject);


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over openpatientdir.
function openpatientdir_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to openpatientdir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ea_gethelp(get(handles.leadfigure,'SelectionType'),hObject);


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over viewmanual.
function viewmanual_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to viewmanual (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ea_gethelp(get(handles.leadfigure,'SelectionType'),hObject);


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over patdir_choosebox.
function patdir_choosebox_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to patdir_choosebox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ea_gethelp(get(handles.leadfigure,'SelectionType'),hObject);


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over recent.
function recent_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to recent (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ea_gethelp(get(handles.leadfigure,'SelectionType'),hObject);


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over electrode_model_popup.
function electrode_model_popup_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to electrode_model_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ea_gethelp(get(handles.leadfigure,'SelectionType'),hObject);


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over MRCT.
function MRCT_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to MRCT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ea_gethelp(get(handles.leadfigure,'SelectionType'),hObject);


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over coregctmethod.
function coregctmethod_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to coregctmethod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ea_gethelp(get(handles.leadfigure,'SelectionType'),hObject);


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over atlassetpopup.
function atlassetpopup_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to atlassetpopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ea_gethelp(get(handles.leadfigure,'SelectionType'),hObject);


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over vizspacepopup.
function vizspacepopup_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to vizspacepopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ea_gethelp(get(handles.leadfigure,'SelectionType'),hObject);


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over specify2dwrite.
function specify2dwrite_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to specify2dwrite (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ea_gethelp(get(handles.leadfigure,'SelectionType'),hObject);


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over cmappushbutton.
function cmappushbutton_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to cmappushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ea_gethelp(get(handles.leadfigure,'SelectionType'),hObject);


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over dlgroupc.
function dlgroupc_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to dlgroupc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ea_gethelp(get(handles.leadfigure,'SelectionType'),hObject);


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over updatebutn.
function updatebutn_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to updatebutn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ea_gethelp(get(handles.leadfigure,'SelectionType'),hObject);


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over coregmrmethod.
function coregmrmethod_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to coregmrmethod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ea_gethelp(get(handles.leadfigure,'SelectionType'),hObject);


% --- Executes on selection change in sidespopup.
function sidespopup_Callback(hObject, eventdata, handles)
% hObject    handle to sidespopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns sidespopup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from sidespopup
ea_gethelp(get(handles.leadfigure,'SelectionType'),hObject);


% --- Executes during object creation, after setting all properties.
function sidespopup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sidespopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in normsettings.
function normsettings_Callback(hObject, eventdata, handles)
% hObject    handle to normsettings (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
normsettingsfunc = getappdata(handles.normsettings,'normsettingsfunc');
feval(normsettingsfunc, handles);


% --- Executes on button press in checkfigures.
function checkfigures_Callback(hObject, eventdata, handles)
% hObject    handle to checkfigures (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkfigures


% --- Executes on button press in scrf.
function scrf_Callback(hObject, eventdata, handles)
% hObject    handle to scrf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of scrf


% --- Executes on selection change in reconmethod.
function reconmethod_Callback(hObject, eventdata, handles)
% hObject    handle to reconmethod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns reconmethod contents as cell array
%        contents{get(hObject,'Value')} returns selected item from reconmethod
if get(hObject,'Value')<3 % TRAC/CORE
    set(handles.targetpopup,'enable','on');
    set(handles.maskwindow_txt,'enable','on');
else % PACER
    prefs=ea_prefs;
    set(handles.targetpopup,'enable','off');
    set(handles.maskwindow_txt,'enable','off');
end

% --- Executes during object creation, after setting all properties.
function reconmethod_CreateFcn(hObject, eventdata, handles)
% hObject    handle to reconmethod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in coreg_checkbox.
function coreg_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to coreg_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of coreg_checkbox

ea_storeui(handles);


% --- Executes on button press in coregmrcheck.
function coregmrcheck_Callback(hObject, eventdata, handles)
% hObject    handle to coregmrcheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of coregmrcheck
ea_storeui(handles);


% --- Executes on button press in overwriteapproved.
function overwriteapproved_Callback(hObject, eventdata, handles)
% hObject    handle to overwriteapproved (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of overwriteapproved
ea_storeui(handles);


% --- Executes on button press in previouspatient.
function previouspatient_Callback(hObject, eventdata, handles)
% hObject    handle to previouspatient (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Select previous patient in same root folder
uipatdir=getappdata(handles.leadfigure,'uipatdir');
if length(uipatdir)>1
   %ea_error('Selecting the previous patient in folder only works if a single patient was selected.');
elseif isempty(uipatdir)
    load([ea_getearoot,'common',filesep,'ea_recentpatients.mat']);
    if iscell(recentfolders)
        recentfolders=recentfolders(1);
    end

    if strcmp('No recent patients found',recentfolders)
        return
    end

    ea_load_pts(handles,recentfolders);
    return
end

[pth,fn]=fileparts(uipatdir{1});

opts=dir(pth);
pts={opts.name};
pts=pts((cell2mat({opts.isdir})));
todel=[];
for pt=1:length(pts)
   if strcmp(pts{pt}(1),'.')
       todel=[todel,pt];
   end
end
pts(todel)=[];
[~,ix]=ismember(fn,pts);
if ix>1
    nuix=ix-1;
else
    nuix=ix;
end

ea_load_pts(handles,{[pth,filesep,pts{nuix}]});
fprintf('\nSwitched to previous patient: %s\n\n', pts{nuix});
if isfield(handles,'atlassetpopup') % not present in connectome mapper
    options.prefs=ea_prefs;
    atlasset=get(handles.atlassetpopup,'String');
    atlasset=atlasset{get(handles.atlassetpopup,'Value')};

    ea_listatlassets(options,handles,get(handles.vizspacepopup,'Value'),atlasset);
end

% --- Executes on button press in nextpatient.
function nextpatient_Callback(hObject, eventdata, handles)
% hObject    handle to nextpatient (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Select next patient in same root folder
uipatdir=getappdata(handles.leadfigure,'uipatdir');
if length(uipatdir)>1 % still works
    %  ea_error('Selecting the next patient in folder only works if a single patient was selected.');
elseif isempty(uipatdir)
    % load recent patient then.

    load([ea_getearoot,'common',filesep,'ea_recentpatients.mat']);
    if iscell(recentfolders)
        recentfolders=recentfolders(1);
    end

    if strcmp('No recent patients found',recentfolders)
        return
    end

    ea_load_pts(handles,recentfolders);
    return
    %   ea_error('Selecting the next patient in folder only works if a patient was selected before.');
end

[pth,fn]=fileparts(uipatdir{1});

opts=dir(pth);
pts={opts.name};
pts=pts((cell2mat({opts.isdir})));
todel=[];
for pt=1:length(pts)
   if strcmp(pts{pt}(1),'.')
       todel=[todel,pt];
   end
end
pts(todel)=[];
[~,ix]=ismember(fn,pts);

if length(pts)>ix
    nuix=ix+1;
else
    nuix=ix;
end
ea_load_pts(handles,{[pth,filesep,pts{nuix}]});
fprintf('\nSwitched to next patient: %s\n\n', pts{nuix});
if isfield(handles,'atlassetpopup') % not present in connectome mapper
    options.prefs=ea_prefs;
    atlasset=get(handles.atlassetpopup,'String');
    atlasset=atlasset{get(handles.atlassetpopup,'Value')};

    ea_listatlassets(options,handles,get(handles.vizspacepopup,'Value'),atlasset);
end

% --- Executes on button press in slicer_original.
function slicer_original_Callback(hObject, eventdata, handles)
% hObject    handle to slicer_original (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
options = ea_handles2options(handles);
options.uipatdirs=getappdata(handles.leadfigure,'uipatdir');
bids = getappdata(handles.leadfigure, 'bids');
subjId = getappdata(handles.leadfigure, 'subjId');
options.subj = bids.getSubj(subjId{1}, options.modality);
ea_runslicer(options, 1);


% --- Executes on button press in slicer_coregistered.
function slicer_coregistered_Callback(hObject, eventdata, handles)
% hObject    handle to slicer_coregistered (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
options = ea_handles2options(handles);
options.uipatdirs=getappdata(handles.leadfigure,'uipatdir');
bids = getappdata(handles.leadfigure, 'bids');
subjId = getappdata(handles.leadfigure, 'subjId');
options.subj = bids.getSubj(subjId{1}, options.modality);
ea_runslicer(options, 2);


% --- Executes on button press in slicer_normalized.
function slicer_normalized_Callback(hObject, eventdata, handles)
% hObject    handle to slicer_normalized (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
options = ea_handles2options(handles);
options.uipatdirs=getappdata(handles.leadfigure,'uipatdir');
bids = getappdata(handles.leadfigure, 'bids');
subjId = getappdata(handles.leadfigure, 'subjId');
options.subj = bids.getSubj(subjId{1}, options.modality);
ea_runslicer(options, 3);


% --- Executes on button press in slicer_contact.
function slicer_contact_Callback(hObject, eventdata, handles)
% hObject    handle to slicer_contact (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
options = ea_handles2options(handles);
options.uipatdirs=getappdata(handles.leadfigure,'uipatdir');
bids = getappdata(handles.leadfigure, 'bids');
subjId = getappdata(handles.leadfigure, 'subjId');
options.subj = bids.getSubj(subjId{1}, options.modality);
ea_runslicer(options, 4);


% --- Executes on button press in side1.
function side1_Callback(hObject, eventdata, handles)
% hObject    handle to side1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of side1


% --- Executes on button press in side2.
function side2_Callback(hObject, eventdata, handles)
% hObject    handle to side2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of side2


% --- Executes on button press in side3.
function side3_Callback(hObject, eventdata, handles)
% hObject    handle to side3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of side3
ea_checktargetentry(handles);

% --- Executes on button press in side4.
function side4_Callback(hObject, eventdata, handles)
% hObject    handle to side4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of side4
ea_checktargetentry(handles);

% --- Executes on button press in side5.
function side5_Callback(hObject, eventdata, handles)
% hObject    handle to side5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of side5
ea_checktargetentry(handles);

% --- Executes on button press in side6.
function side6_Callback(hObject, eventdata, handles)
% hObject    handle to side6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of side6
ea_checktargetentry(handles);

% --- Executes on button press in side7.
function side7_Callback(hObject, eventdata, handles)
% hObject    handle to side7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of side7
ea_checktargetentry(handles);

% --- Executes on button press in side8.
function side8_Callback(hObject, eventdata, handles)
% hObject    handle to side8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of side8
ea_checktargetentry(handles);

% --- Executes on button press in side9.
function side9_Callback(hObject, eventdata, handles)
% hObject    handle to side9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of side9
ea_checktargetentry(handles);

% --- Executes on button press in side10.
function side10_Callback(hObject, eventdata, handles)
% hObject    handle to side10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of side10
ea_checktargetentry(handles);


% --- Executes on button press in side11.
function side11_Callback(hObject, eventdata, handles)
% hObject    handle to side11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of side11
ea_checktargetentry(handles);


% --- Executes on button press in side12.
function side12_Callback(hObject, eventdata, handles)
% hObject    handle to side12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of side12
ea_checktargetentry(handles);


% --- Executes on button press in side13.
function side13_Callback(hObject, eventdata, handles)
% hObject    handle to side13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of side13
ea_checktargetentry(handles);


% --- Executes on button press in side14.
function side14_Callback(hObject, eventdata, handles)
% hObject    handle to side14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of side14
ea_checktargetentry(handles);


% --- Executes on button press in side15.
function side15_Callback(hObject, eventdata, handles)
% hObject    handle to side15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of side15
ea_checktargetentry(handles);


function ea_checktargetentry(handles)
higherelchosen=0;
elnum = sum(cellfun(@(f) ~isempty(f), regexp(fieldnames(handles),'^side\d+$','match')));
for el=3:elnum
    if get(handles.(['side',num2str(el)]),'Value')
        higherelchosen=1;
    end
end
if higherelchosen
    set(handles.targetpopup,'Value',3);
    set(handles.targetpopup,'Enable','off');
else
    if get(handles.reconmethod,'Value')==1
        set(handles.targetpopup,'Enable','on');
    end
end


% --- Executes on selection change in scrfmask.
function scrfmask_Callback(hObject, eventdata, handles)
% hObject    handle to scrfmask (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns scrfmask contents as cell array
%        contents{get(hObject,'Value')} returns selected item from scrfmask


% --- Executes during object creation, after setting all properties.
function scrfmask_CreateFcn(hObject, eventdata, handles)
% hObject    handle to scrfmask (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in dcm2niiselect.
function dcm2niiselect_Callback(hObject, eventdata, handles)
% hObject    handle to dcm2niiselect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns dcm2niiselect contents as cell array
%        contents{get(hObject,'Value')} returns selected item from dcm2niiselect


% --- Executes during object creation, after setting all properties.
function dcm2niiselect_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dcm2niiselect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in extractsurface.
function extractsurface_Callback(hObject, eventdata, handles)
% hObject    handle to extractsurface (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of extractsurface


% --- Executes on button press in localizeecog.
function localizeecog_Callback(hObject, eventdata, handles)
% hObject    handle to localizeecog (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of localizeecog


% --- Executes on selection change in surfacemethod.
function surfacemethod_Callback(hObject, eventdata, handles)
% hObject    handle to surfacemethod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns surfacemethod contents as cell array
%        contents{get(hObject,'Value')} returns selected item from surfacemethod
if strcmpi(hObject.String{hObject.Value}, 'FreeSurfer')
    handles.surfsettings.Visible = 'on';
else
    handles.surfsettings.Visible = 'off';
end


% --- Executes during object creation, after setting all properties.
function surfacemethod_CreateFcn(hObject, eventdata, handles)
% hObject    handle to surfacemethod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in checkcoreg.
function checkreg_Callback(hObject, eventdata, handles)
% hObject    handle to checkcoreg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkcoreg

% --- Executes on button press in refinefit.
function refinefit_Callback(hObject, eventdata, handles)
% hObject    handle to refinefit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of refinefit


% --- Executes on button press in surfsettings.
function surfsettings_Callback(hObject, eventdata, handles)
% hObject    handle to surfsettings (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ea_freesurfersetting;
