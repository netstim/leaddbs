function varargout = lead_anatomy(varargin)
% LEAD_ANATOMY MATLAB code for lead_anatomy.fig
%      LEAD_ANATOMY, by itself, creates a new LEAD_ANATOMY or raises the existing
%      singleton*.
%
%      H = LEAD_ANATOMY returns the handle to a new LEAD_ANATOMY or the handle to
%      the existing singleton*.
%
%      LEAD_ANATOMY('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in LEAD_ANATOMY.M with the given input arguments.
%
%      LEAD_ANATOMY('Property','Value',...) creates a new LEAD_ANATOMY or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before lead_anatomy_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to lead_anatomy_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help lead_anatomy

% Last Modified by GUIDE v2.5 02-Mar-2023 19:13:49

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @lead_anatomy_OpeningFcn, ...
                   'gui_OutputFcn',  @lead_anatomy_OutputFcn, ...
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


% --- Executes just before lead_anatomy is made visible.
function lead_anatomy_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to lead_anatomy (see VARARGIN)

handles.prod = 'anatomy';
handles.callingfunction = 'lead_anatomy';

% add recentpatients patients...
ea_initrecent(handles, 'patients');
earoot=ea_getearoot;
if ~isdeployed
    addpath(genpath(earoot));
    rmpath(genpath([earoot,'.git']));
    rmpath(genpath([earoot,'release']));
end

% for now disable space dropdown
set(handles.vizspacepopup,'enable','off');

options.prefs = ea_prefs('');

% load atlassets
ea_listatlassets(options, handles, 1);

% list templates
list=ea_assignbackdrop('list',options,'Patient');
set(handles.tdbackdrop,'String',list);

set(hObject,'Color',[1 1 1]);
set(handles.versiontxt,'String',['v',ea_getvsn('local')]);

% add logo
set(0,'CurrentFigure',handles.leadfigure);
im=imread([earoot,'icons',filesep,'logo_lead_anatomy.png']);
image(im);
axis off;
axis equal;
set(handles.leadfigure,'name','Lead Anatomy','color','w');

ea_init_coregmrpopup(handles, options.prefs.mrcoreg.default);

% Initialize norm methods popupmenu
ea_init_normpopup(handles, options.prefs.normalize.default);

ea_processguiargs(handles,varargin)

% add tools menu
ea_menu_initmenu(handles,{'export','cluster','prefs','transfer','space'},options.prefs);

ea_firstrun(handles,options);

ea_bind_dragndrop(handles.leadfigure, ...
    @(obj,evt) DropFcn(obj,evt,handles), ...
    @(obj,evt) DropFcn(obj,evt,handles));

% Choose default command line output for lead_anatomy
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes lead_anatomy wait for user response (see UIRESUME)
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

ea_busyaction('on',handles.leadfigure,'anatomy');
if ~isempty(patdir)
    ea_load_pts(handles, patdir);

    % refresh backdrop list if only one patient dragged in (use standard list for mutiple patients)
    if length(patdir) == 1
        options.prefs=ea_prefs('');
        [options.root,options.patientname]=fileparts(get(handles.patdir_choosebox,'String'));
        options.root=[options.root,filesep];
        list=ea_assignbackdrop('list',options,'Patient');
        set(handles.tdbackdrop,'String',list);
    end

    % disable atlas list refreshing for now
%     if isfield(handles,'atlassetpopup')
%         atlasset=get(handles.atlassetpopup,'String');
%         atlasset=atlasset{get(handles.atlassetpopup,'Value')};
%         options.prefs=ea_prefs('');
%         ea_listatlassets(options,handles,get(handles.vizspacepopup,'Value'),atlasset);
%     end
end
ea_busyaction('off',handles.leadfigure,'anatomy');


% --- Outputs from this function are returned to the command line.
function varargout = lead_anatomy_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in viz3d.
function viz3d_Callback(hObject, eventdata, handles)
% hObject    handle to viz3d (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

leadfig=handles.leadfigure;
ea_busyaction('on',leadfig,'anatomy');

options=ea_handles2options(handles);
options.macaquemodus=0;

options.uipatdirs=getappdata(handles.leadfigure,'uipatdir');

options.d3.write=1;

options.leadprod = 'anatomy';

ea_run('run',options);

ea_busyaction('off',leadfig,'anatomy');


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


% --- Executes on selection change in vizspacepopup.
function vizspacepopup_Callback(hObject, eventdata, handles)
% hObject    handle to vizspacepopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns vizspacepopup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from vizspacepopup


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


% --- Executes on button press in patdir_choosebox.
function patdir_choosebox_Callback(hObject, eventdata, handles)
% hObject    handle to patdir_choosebox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
options.prefs=ea_prefs('');
ea_getpatients(options,handles);
[options.root,options.patientname]=fileparts(get(handles.patdir_choosebox,'String'));
options.root=[options.root,filesep];
list=ea_assignbackdrop('list',options,'Patient');
set(handles.tdbackdrop,'String',list);


% --- Executes on selection change in recentpatients.
function recentpatients_Callback(hObject, eventdata, handles)
% hObject    handle to recentpatients (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns recentpatients contents as cell array
%        contents{get(hObject,'Value')} returns selected item from recentpatients
ea_recentcallback(handles, 'patients');
options.prefs=ea_prefs('');
[options.root,options.patientname]=fileparts(get(handles.patdir_choosebox,'String'));
options.root=[options.root,filesep];
list=ea_assignbackdrop('list',options,'Patient');
set(handles.tdbackdrop,'String',list);


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


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in tracor.
function tracor_Callback(hObject, eventdata, handles)
% hObject    handle to tracor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns tracor contents as cell array
%        contents{get(hObject,'Value')} returns selected item from tracor
switch get(hObject,'Value')
    case 1
        set(handles.xyzslice,'String','z = ');
    case 2
        set(handles.xyzslice,'String','y = ');
    case 3
        set(handles.xyzslice,'String','x = ');
end

% --- Executes during object creation, after setting all properties.
function tracor_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tracor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function depth_Callback(hObject, eventdata, handles)
% hObject    handle to xdepth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of xdepth as text
%        str2double(get(hObject,'String')) returns contents of xdepth as a double


% --- Executes during object creation, after setting all properties.
function depth_CreateFcn(hObject, eventdata, handles)
% hObject    handle to xdepth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in viz2d.
function viz2d_Callback(hObject, eventdata, handles)
% hObject    handle to viz2d (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
leadfig=handles.leadfigure;
ea_busyaction('on',leadfig,'anatomy');

options=ea_handles2options(handles);
options.macaquemodus=0;

options.d2=ea_tdhandles2options(handles);

% ea_setprefs('d2',d2);

options.uipatdirs=getappdata(handles.leadfigure,'uipatdir');
if isempty(options.uipatdirs)
    options.uipatdirs={''};
end

if strcmp(options.d2.backdrop, 'Choose...')
    options.d2.backdrop = getappdata(gcf,'customfile');
end

for pt=1:length(options.uipatdirs)
    [pth,ptname]=fileparts(options.uipatdirs{pt});
    options.root=[pth,filesep];
    options.patientname=ptname;
    options.prefs=ea_prefs(options.patientname);
%     options.d2=options.prefs.machine.d2;
    options.d2.writeatlases=1;
    options.d2.atlasopacity=0.2;
    options.d2.tracor=get(handles.tracor,'Value');
    options.d2.depth=[str2double(get(handles.depth,'String')),...
        str2double(get(handles.depth,'String')),...
        str2double(get(handles.depth,'String'))];
    setzero=[1,1,1]; setzero(ea_getdims(options.d2.tracor,1))=0;
    options.d2.depth(logical(setzero))=0;
    options.d2.showlegend=0;
    [Vtra,Vcor,Vsag]=ea_assignbackdrop(options.d2.backdrop,options,'Patient');
    Vs={Vtra,Vcor,Vsag};
    options.sides=1;
    h=ea_writeplanes(options,options.d2.depth,options.d2.tracor,Vs{options.d2.tracor},'on',2);
    set(h,'Position',[0,0,800,800]);
end

ea_busyaction('off',leadfig,'anatomy');


% --- Executes on selection change in tdbackdrop.
function tdbackdrop_Callback(hObject, eventdata, handles)
% hObject    handle to tdbackdrop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns tdbackdrop contents as cell array
%        contents{get(hObject,'Value')} returns selected item from tdbackdrop
popvals=get(hObject,'String');
if strcmp(popvals{get(hObject,'Value')},'Choose...')
    [FileName,PathName] = uigetfile('*.nii','Choose anatomical image...');
    if all([PathName,FileName])
        setappdata(gcf,'customfile',[PathName,FileName]);
    else
        set(hObject,'Value',1);
    end
end


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


% --- Executes on button press in tdlabelcheck.
function tdlabelcheck_Callback(hObject, eventdata, handles)
% hObject    handle to tdlabelcheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tdlabelcheck


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


% --- Executes when user attempts to close leadfigure.
function leadfigure_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to leadfigure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure

% d2=ea_tdhandles2options(handles);
% ea_setprefs('d2',d2);

delete(hObject);


% --- Executes on button press in updatebutn.
function updatebutn_Callback(hObject, eventdata, handles)
% hObject    handle to updatebutn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in MRCT.
function MRCT_Callback(hObject, eventdata, handles)
% hObject    handle to MRCT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns MRCT contents as cell array
%        contents{get(hObject,'Value')} returns selected item from MRCT
ea_switchctmr(handles,get(hObject,'Value'));


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


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over patdir_choosebox.
function patdir_choosebox_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to patdir_choosebox (see GCBO)
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


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over coregmrmethod.
function coregmrmethod_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to coregmrmethod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ea_gethelp(get(handles.leadfigure,'SelectionType'),hObject);


% --- Executes on button press in tdfidcheck.
function tdfidcheck_Callback(hObject, eventdata, handles)
% hObject    handle to tdfidcheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tdfidcheck


% --- Executes on button press in specify2dwrite.
function specify2dwrite_Callback(hObject, eventdata, handles)
% hObject    handle to specify2dwrite (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
options.native=get(handles.vizspacepopup,'Value')==2;
[options.root,options.patientname]=fileparts(get(handles.patdir_choosebox,'String'));
options.root=[options.root,filesep];
options.modality=get(handles.MRCT,'Value');
options.prefs=ea_prefs(options.patientname);
ea_spec2dwrite(options);
