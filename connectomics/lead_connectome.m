function varargout = lead_connectome(varargin)
% LEADFIGURE MATLAB code for leadfigure.fig
%      LEADFIGURE, by itself, creates a new LEADFIGURE or raises the existing
%      singleton*.
%
%      H = LEADFIGURE returns the handle to a new LEADFIGURE or the handle to
%      the existing singleton*.
%
%      LEADFIGURE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in LEADFIGURE.M with the given input arguments.
%
%      LEADFIGURE('Property','Value',...) creates a new LEADFIGURE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before lead_connectome_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to lead_connectome_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help leadfigure

% Last Modified by GUIDE v2.5 02-Mar-2023 19:11:26

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @lead_connectome_OpeningFcn, ...
                   'gui_OutputFcn',  @lead_connectome_OutputFcn, ...
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


% --- Executes just before leadfigure is made visible.
function lead_connectome_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to leadfigure (see VARARGIN)

earoot=ea_getearoot;
im=imread([earoot,'icons',filesep,'logo_lead_connectome.png']);
image(im);
axes(handles.logoaxes);
axis off;
axis equal;
set(handles.leadfigure,'name','Lead Connectome','color','w');

% Add parcellations to menu:
parcFiles = dir([ea_space([],'labeling'),'*.nii']);
if ~isempty(parcFiles)
    parcellations = cellfun(@(x) {strrep(x, '.nii', '')}, {parcFiles.name});
    set(handles.parcellation,'String',parcellations);

    options.prefs = ea_prefs;
    parc = find(ismember(parcellations,options.prefs.lc.defaultParcellation));
    if ~isempty(parc)
        set(handles.parcellation,'Value',parc);
    end
end

set(handles.versiontxt,'String',['v',ea_getvsn('local')]);

% add ft methods to menu
cnt=1;
ndir=dir([earoot,'connectomics',filesep,'ea_ft_*.m']);
ftFunctions=cell(0);
ftMethods=cell(0);
for nd=length(ndir):-1:1
    [~,methodf]=fileparts(ndir(nd).name);
    try
        [thisndc,spmvers]=eval([methodf,'(','''prompt''',')']);
        if ismember(spm('ver'),spmvers)
            ftMethods{cnt}=thisndc;
            ftFunctions{cnt}=methodf;
            cnt=cnt+1;
        end
    end
end

setappdata(gcf,'ftFunctions',ftFunctions);
set(handles.ftmethod,'String',ftMethods);

% Initialize norm methods popupmenu
ea_init_normpopup(handles, options.prefs.normalize.default);

% add recentpatients patients...
ea_initrecent(handles, 'patients');

% update UI:
try
    lc=options.prefs.machine.lc;
catch
    lc=ea_initlcopts(handles.leadfigure);
end
lc2handles(lc,handles);

ea_init_coregmrpopup(handles, options.prefs.mrcoreg.default);

if isempty(varargin) % "standard alone" mode, i.e. not dependend from lead
    isindependent=1;
else
    if strcmp(varargin{1},'dependent')
        isindependent=0;
        ea_makehidelcindep(handles);
    else
        isindependent=1;
    end
end
setappdata(handles.leadfigure,'isindependent',isindependent);

if isindependent
    handles.prod='connectome';
else
    handles.prod='dbs_connectome';
end
handles.callingfunction='lead_connectome';


if isindependent
    %% add tools menu
    ea_processguiargs(handles,varargin)

    ea_menu_initmenu(handles,{'export','cluster','prefs','transfer','space','methods'},ea_prefs);

end

ea_bind_dragndrop(handles.leadfigure, ...
    @(obj,evt) DropFcn(obj,evt,handles), ...
    @(obj,evt) DropFcn(obj,evt,handles));

ea_ListBoxRenderer(handles.parcellation);

ea_firstrun(handles,options);

% Choose default command line output for leadfigure
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);


% UIWAIT makes leadfigure wait for user response (see UIRESUME)
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

ea_busyaction('on',handles.leadfigure,'connectome');
if ~isempty(patdir)
    ea_load_pts(handles, patdir);

    if isfield(handles,'atlassetpopup')
        atlasset=get(handles.atlassetpopup,'String');
        atlasset=atlasset{get(handles.atlassetpopup,'Value')};
        options.prefs=ea_prefs('');
        ea_listatlassets(options,handles,get(handles.vizspacepopup,'Value'),atlasset);
    end
end
ea_busyaction('off',handles.leadfigure,'connectome');


function ea_makehidelcindep(handles)

set(handles.importpanel,'visible','off');
set(handles.runsavebutn,'String','Save and close');
set(handles.runsavebutn,'Position',[218,8,154,41]);
set(handles.exportcode,'visible','off');
set(handles.overwriteapproved,'visible','off');
set(handles.openpatientdir,'visible','off');
set(handles.previouspatient,'visible','off');
set(handles.nextpatient,'visible','off');

headermovedown = 238;
panelmovedown = 22;

pos=get(handles.logoaxes,'Position');
pos(2)=pos(2)-headermovedown;
set(handles.logoaxes,'Position',pos);

pos=get(handles.titletext,'Position');
pos(2)=pos(2)-headermovedown;
set(handles.titletext,'Position',pos);

pos=get(handles.versiontxt,'Position');
pos(2)=pos(2)-headermovedown;
set(handles.versiontxt,'Position',pos);

pos=get(handles.nclinuse,'Position');
pos(2)=pos(2)-headermovedown;
set(handles.nclinuse,'Position',pos);

pos=get(handles.genpanel,'Position');
pos(2)=pos(2)-panelmovedown;
set(handles.genpanel,'Position',pos);

pos=get(handles.strucpanel,'Position');
pos(2)=pos(2)-panelmovedown;
set(handles.strucpanel,'Position',pos);

pos=get(handles.funcpanel,'Position');
pos(2)=pos(2)-panelmovedown;
set(handles.funcpanel,'Position',pos);

pos=get(handles.graphpanel,'Position');
pos(2)=pos(2)-panelmovedown;
set(handles.graphpanel,'Position',pos);

pos=get(handles.leadfigure,'Position');
pos(4)=pos(4)-headermovedown;
set(handles.leadfigure,'Position',pos);


% --- Outputs from this function are returned to the command line.
function varargout = lead_connectome_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1}=hObject;

% --- Executes on selection change in parcellation.
function parcellation_Callback(hObject, eventdata, handles)
% hObject    handle to parcellation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns parcellation contents as cell array
%        contents{get(hObject,'Value')} returns selected item from parcellation


% --- Executes during object creation, after setting all properties.
function parcellation_CreateFcn(hObject, eventdata, handles)
% hObject    handle to parcellation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in compute_CM_func.
function compute_CM_func_Callback(hObject, eventdata, handles)
% hObject    handle to compute_CM_func (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of compute_CM_func


% --- Executes on button press in perf_ft.
function perf_ft_Callback(hObject, eventdata, handles)
% hObject    handle to perf_ft (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of perf_ft


% --- Executes on selection change in ftmethod.
function ftmethod_Callback(hObject, eventdata, handles)
% hObject    handle to ftmethod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns ftmethod contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ftmethod

ftmethod=get(handles.ftmethod,'String');
ftmethod=ftmethod{get(handles.ftmethod,'Value')};
[gqiUItxt]=eval(['ea_ft_gqi_yeh','(','''prompt''',')']);

if strcmp(ftmethod, gqiUItxt)
    set(handles.fiber_count, 'Visible', 'on');
    set(handles.fiber_count_txt, 'Visible', 'on');

    if ismember(get(handles.upsamplingfactor,'value'),[3,5])
        set(handles.use_internal_upsampling,'enable','on');
    else
        set(handles.use_internal_upsampling,'enable','off');
        set(handles.use_internal_upsampling,'Value',0);
    end

else
    set(handles.use_internal_upsampling,'enable','off');
    set(handles.use_internal_upsampling,'Value',0);
    set(handles.fiber_count, 'Visible', 'off');
    set(handles.fiber_count_txt, 'Visible', 'off');
end



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


% --- Executes on button press in compute_CM_struc.
function compute_CM_struc_Callback(hObject, eventdata, handles)
% hObject    handle to compute_CM_struc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of compute_CM_struc


% --- Executes on button press in runsavebutn.
function runsavebutn_Callback(hObject, eventdata, handles)
% hObject    handle to runsavebutn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

ea_busyaction('on', handles.leadfigure, 'connectome');

options = ea_handles2options(handles);
options.uipatdirs = getappdata(handles.leadfigure,'uipatdir');

isindependent = getappdata(handles.leadfigure,'isindependent');
ea_savelcopts(handles)

% run execution:
if isindependent
    options.leadprod = 'connectome';
    ea_run('run',options);
end

ea_busyaction('off', handles.leadfigure, 'connectome');


% --- Executes on button press in degree_centrality.
function degree_centrality_Callback(hObject, eventdata, handles)
% hObject    handle to degree_centrality (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of degree_centrality


% --- Executes on button press in eigenvector_centrality.
function eigenvector_centrality_Callback(hObject, eventdata, handles)
% hObject    handle to eigenvector_centrality (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of eigenvector_centrality


% --- Executes on button press in nodal_efficiency.
function nodal_efficiency_Callback(hObject, eventdata, handles)
% hObject    handle to nodal_efficiency (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of nodal_efficiency


% --- Executes on button press in struc_func_sim.
function struc_func_sim_Callback(hObject, eventdata, handles)
% hObject    handle to struc_func_sim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of struc_func_sim


% --- Executes on button press in compute_GM_func.
function compute_GM_func_Callback(hObject, eventdata, handles)
% hObject    handle to compute_GM_func (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of compute_GM_func


% --- Executes on button press in compute_GM_struc.
function compute_GM_struc_Callback(hObject, eventdata, handles)
% hObject    handle to compute_GM_struc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of compute_GM_struc


function handles=lc2handles(lc,handles)

% General settings
parcellations = get(handles.parcellation,'String');
if ~isfield(lc.general, 'parcellation')
    options.prefs = ea_prefs;
    defaultParc = options.prefs.lc.defaultParcellation;
    set(handles.parcellation,'Value',find(ismember(parcellations, defaultParc)));
else
    parcIdx = find(ismember(parcellations, lc.general.parcellation), 1);
    if ~isempty(parcIdx)
        set(handles.parcellation,'Value',parcIdx);
    else
        options.prefs = ea_prefs;
        parc = find(ismember(parcellations, options.prefs.lc.defaultParcellation));
        if ~isempty(parc)
            set(handles.parcellation,'Value',parc);
        end
    end
end

% Graph options:
try set(handles.struc_func_sim,'Value',lc.graph.struc_func_sim); end
try set(handles.nodal_efficiency,'Value',lc.graph.nodal_efficiency); end
try set(handles.eigenvector_centrality,'Value',lc.graph.eigenvector_centrality); end
try set(handles.degree_centrality,'Value',lc.graph.degree_centrality); end
try set(handles.fthresh,'String',num2str(lc.graph.fthresh)); end
try set(handles.sthresh,'String',num2str(lc.graph.sthresh)); end


% functional options:
try set(handles.compute_CM_func,'Value',lc.func.compute_CM); end
try set(handles.compute_GM_func,'Value',lc.func.compute_GM); end
try set(handles.TR,'String',num2str(lc.func.prefs.TR)); end


% structural options:
try set(handles.compute_CM_struc,'Value',lc.struc.compute_CM); end
try set(handles.compute_GM_struc,'Value',lc.struc.compute_GM); end
ftFunctions = getappdata(handles.leadfigure, 'ftFunctions');
ftMethodIdx = find(ismember(ftFunctions, lc.struc.ft.method));
if ~isempty(ftMethodIdx)
    set(handles.ftmethod,'Value',ftMethodIdx);
else
    options.prefs = ea_prefs;
    defaultftMethod = options.prefs.machine.lc.struc.ft.method;
    set(handles.ftmethod,'Value',find(ismember(ftFunctions, defaultftMethod)));
end

if strcmp(lc.struc.ft.method, 'ea_ft_gqi_yeh')
    try set(handles.fiber_count, 'Visible', 'on'); end
    try set(handles.fiber_count_txt, 'Visible', 'on'); end
else
    try set(handles.fiber_count, 'Visible', 'off'); end
    try set(handles.fiber_count_txt, 'Visible', 'off'); end
end

try set(handles.fiber_count, 'String', num2str(lc.struc.ft.dsistudio.fiber_count)); end

try set(handles.normalize_fibers,'Value',lc.struc.ft.normalize); end
try set(handles.perf_ft,'Value',lc.struc.ft.do); end


% --- Executes on button press in normalize_fibers.
function normalize_fibers_Callback(hObject, eventdata, handles)
% hObject    handle to normalize_fibers (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of normalize_fibers


function TR_Callback(hObject, eventdata, handles)
% hObject    handle to TR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of TR as text
%        str2double(get(hObject,'String')) returns contents of TR as a double


% --- Executes during object creation, after setting all properties.
function TR_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double


% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function sthresh_Callback(hObject, eventdata, handles)
% hObject    handle to sthresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sthresh as text
%        str2double(get(hObject,'String')) returns contents of sthresh as a double


% --- Executes during object creation, after setting all properties.
function sthresh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sthresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function fthresh_Callback(hObject, eventdata, handles)
% hObject    handle to fthresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fthresh as text
%        str2double(get(hObject,'String')) returns contents of fthresh as a double


% --- Executes during object creation, after setting all properties.
function fthresh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fthresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in normalize_checkbox.
function normalize_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to normalize_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of normalize_checkbox


% --- Executes on selection change in normmethod.
function normmethod_Callback(hObject, eventdata, handles)
% hObject    handle to normmethod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns normmethod contents as cell array
%        contents{get(hObject,'Value')} returns selected item from normmethod
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


% --- Executes on button press in normcheck.
function normcheck_Callback(hObject, eventdata, handles)
% hObject    handle to normcheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of normcheck


% --- Executes on button press in exportcode.
function exportcode_Callback(hObject, eventdata, handles)
% hObject    handle to exportcode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

ea_busyaction('on', handles.leadfigure, 'connectome');

ea_savelcopts(handles);

options = ea_handles2options(handles);
options.uipatdirs = getappdata(handles.leadfigure,'uipatdir');
options.macaquemodus = 0;

options.leadprod = 'connectome';

ea_run('export',options);

ea_busyaction('off', handles.leadfigure, 'connectome');


% --- Executes on button press in patdir_choosebox.
function patdir_choosebox_Callback(hObject, eventdata, handles)
% hObject    handle to patdir_choosebox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ea_busyaction('on',handles.leadfigure,'connectome');
options.prefs=ea_prefs('');
ea_getpatients(options,handles);
ea_busyaction('off',handles.leadfigure,'connectome');

% --- Executes on selection change in recentpatients.
function recentpatients_Callback(hObject, eventdata, handles)
% hObject    handle to recentpatients (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns recentpatients contents as cell array
%        contents{get(hObject,'Value')} returns selected item from recentpatients
ea_busyaction('on',handles.leadfigure,'connectome');
ea_recentcallback(handles, 'patients');
ea_busyaction('off',handles.leadfigure,'connectome');


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


% --- Executes on button press in include_lead_connectome_subroutine.
function include_lead_connectome_subroutine_Callback(hObject, eventdata, handles)
% hObject    handle to include_lead_connectome_subroutine (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of include_lead_connectome_subroutine


% --- Executes on button press in coreg_checkbox.
function coreg_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to coreg_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of coreg_checkbox
ea_storeui(handles);


% --- Executes on selection change in coregmrmethod.
function coregmrmethod_Callback(hObject, eventdata, handles)
% hObject    handle to coregmrmethod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns coregmrmethod contents as cell array
%        contents{get(hObject,'Value')} returns selected item from coregmrmethod
prefs = ea_prefs;
if strcmp(hObject.String{hObject.Value}, 'ANTs') &&  prefs.env.dev
    set(handles.addSyN, 'Visible', 'on');
else
    set(handles.addSyN, 'Visible', 'off');
end


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


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over coregmrmethod.
function coregmrmethod_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to coregmrmethod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ea_gethelp(get(handles.leadfigure,'SelectionType'),hObject);


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


% --- Executes on button press in openpatientdir.
function openpatientdir_Callback(hObject, eventdata, handles)
% hObject    handle to openpatientdir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ea_openpatdir(handles);


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
    % load recentpatients patient then.

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
if isfield(handles,'atlassetpopup') % not present in connectome mapper
    options.prefs=ea_prefs;
    atlasset=get(handles.atlassetpopup,'String');
    atlasset=atlasset{get(handles.atlassetpopup,'Value')};

    ea_listatlassets(options,handles,get(handles.vizspacepopup,'Value'),atlasset);
end


% --- Executes on button press in checkregfmri.
function checkregfmri_Callback(hObject, eventdata, handles)
% hObject    handle to checkregfmri (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkregfmri


% --- Executes on button press in checkregdmri.
function checkregdmri_Callback(hObject, eventdata, handles)
% hObject    handle to checkregdmri (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkregdmri


function fiber_count_Callback(hObject, eventdata, handles)
% hObject    handle to fiber_count (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fiber_count as text
%        str2double(get(hObject,'String')) returns contents of fiber_count as a double


% --- Executes during object creation, after setting all properties.
function fiber_count_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fiber_count (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when user attempts to close leadfigure.
function leadfigure_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to leadfigure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure

ea_savelcopts(handles);
delete(hObject);


% --- Executes on selection change in upsamplingfactor.
function upsamplingfactor_Callback(hObject, eventdata, handles)
% hObject    handle to upsamplingfactor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns upsamplingfactor contents as cell array
%        contents{get(hObject,'Value')} returns selected item from upsamplingfactor

ftmethod=get(handles.ftmethod,'String');
ftmethod=ftmethod{get(handles.ftmethod,'Value')};
[gqiUItxt]=eval(['ea_ft_gqi_yeh','(','''prompt''',')']);

if ismember(get(hObject,'value'),[3,5]) && strcmp(gqiUItxt,ftmethod)
    set(handles.use_internal_upsampling,'enable','on');
else
    set(handles.use_internal_upsampling,'enable','off');
    set(handles.use_internal_upsampling,'Value',0);
end


% --- Executes during object creation, after setting all properties.
function upsamplingfactor_CreateFcn(hObject, eventdata, handles)
% hObject    handle to upsamplingfactor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in use_internal_upsampling.
function use_internal_upsampling_Callback(hObject, eventdata, handles)
% hObject    handle to use_internal_upsampling (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of use_internal_upsampling


% --- Executes on button press in addSyN.
function addSyN_Callback(hObject, eventdata, handles)
% hObject    handle to addSyN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of addSyN
