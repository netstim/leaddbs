function varargout = lead_group_connectome(varargin)
% LEAD_GROUP_CONNECTOME MATLAB code for lead_group_connectome.fig
%      LEAD_GROUP_CONNECTOME, by itself, creates a new LEAD_GROUP_CONNECTOME or raises the existing
%      singleton*.
%
%      H = LEAD_GROUP_CONNECTOME returns the handle to a new LEAD_GROUP_CONNECTOME or the handle to
%      the existing singleton*.
%
%      LEAD_GROUP_CONNECTOME('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in LEAD_GROUP_CONNECTOME.M with the given input arguments.
%
%      LEAD_GROUP_CONNECTOME('Property','Value',...) creates a new LEAD_GROUP_CONNECTOME or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before lead_group_connectome_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to lead_group_connectome_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help lead_group_connectome

% Last Modified by GUIDE v2.5 29-Jan-2020 09:10:42

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @lead_group_connectome_OpeningFcn, ...
    'gui_OutputFcn',  @lead_group_connectome_OutputFcn, ...
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


% --- Executes just before lead_group_connectome is made visible.
function lead_group_connectome_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to lead_group_connectome (see VARARGIN)

% Choose default command line output for lead_group_connectome
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes lead_group_connectome wait for user response (see UIRESUME)
% uiwait(handles.leadfigure);

options.earoot = ea_getearoot;
options.prefs = ea_prefs;
setappdata(handles.leadfigure,'earoot',options.earoot);

% set background image
set(gcf,'color','w');
im=imread([options.earoot,'icons',filesep,'logo_lead_group.png']);
image(im);
axis off;
axis equal;

try
    priorselection=find(ismember(fiberscell,stimparams.usefiberset)); % retrieve prior selection of fiberset.
    set(handles.fiberspopup,'Value',priorselection);
catch    % reinitialize using third entry.
    set(handles.fiberspopup,'Value',1);
end

if get(handles.fiberspopup,'Value')>length(get(handles.fiberspopup,'String'))
    set(handles.fiberspopup,'Value',length(get(handles.fiberspopup,'String')));
end

% Labels:
labeling = dir([ea_space(options,'labeling'),'*.nii']);
labeling = cellfun(@(x) {strrep(x, '.nii', '')}, {labeling.name});

set(handles.labelpopup,'String', labeling);

% Initialize parcellation popupmenu
defaultParc = options.prefs.lg.defaultParcellation;
set(handles.labelpopup,'Value',find(ismember(labeling, defaultParc)));

% Set connectome popup
modlist = ea_genmodlist([],[],options,'dmri');
modlist{end+1}='Do not calculate connectivity stats';
set(handles.fiberspopup,'String',modlist);
set(handles.fiberspopup,'Value',length(modlist));

% set version text:
set(handles.versiontxt,'String',['v',ea_getvsn('local')]);

% make listboxes multiselectable:
set(handles.patientlist,'Max',100,'Min',0);
set(handles.grouplist,'Max',100,'Min',0);
set(handles.clinicallist,'Max',100,'Min',0);

M=getappdata(gcf,'M');
if isempty(M)
    % initialize variable M
    M=ea_initializeM_connectome;
end
setappdata(gcf,'M',M);
ea_refresh_lg_connectome(handles);

handles.prod='group';
handles.callingfunction='lead_group_connectome';

ea_firstrun(handles,options);

ea_menu_initmenu(handles,{'prefs','transfer'},options.prefs);

ea_processguiargs(handles,varargin)

ea_bind_dragndrop(handles.leadfigure, ...
    @(obj,evt) DropFcn(obj,evt,handles), ...
    @(obj,evt) DropFcn(obj,evt,handles));


% --- Drag and drop callback to load patdirs.
function DropFcn(~, event, handles)

% check if dropping area is in patient listbox
if event.Location.getX < 325 && event.Location.getX > 24 && ...
   event.Location.getY < 322 && event.Location.getY > 137
    target = 'patientList';
else
    target = 'groupDir';
end

switch event.DropType
    case 'file'
        folders = event.Data;
    case 'string'
        folders = {event.Data};
end

if strcmp(target, 'groupDir')
    if length(folders) > 1 || ~exist(folders{1}, 'dir')
        ea_error('To choose the group analysis directory, please drag a single folder into Lead Group!', simpleStack = 1);
    end

    groupdir = [folders{1}, filesep];
    set(handles.groupdir_choosebox, 'String', groupdir);
    set(handles.groupdir_choosebox, 'TooltipString', groupdir);

    ea_busyaction('on',handles.leadfigure,'group');

    M = ea_initializeM_connectome;
    M.root = groupdir;

    try % if file already exists, load it (and overwrite M).
        load([groupdir,'GroupConnectomeAnalysis.mat']);
    catch % if not, store it saving M.
        save([groupdir,'GroupConnectomeAnalysis.mat'],'M','-v7.3');
    end

    setappdata(handles.leadfigure,'M',M);

    ea_busyaction('off',handles.leadfigure,'group');

    ea_refresh_lg_connectome(handles);
else
    if strcmp(handles.groupdir_choosebox.String, 'Choose Group Directory')
        ea_error('Please choose a group directory first to store the group analysis!', simpleStack = 1);
    end

    nonexist = cellfun(@(x) ~exist(x, 'dir'), folders);
    if any(nonexist)
        fprintf('\nExcluded non-existent/invalid folder:\n');
        cellfun(@disp, folders(nonexist));
        fprintf('\n');
        folders(nonexist) = [];
    end

    if ~isempty(folders)
        M=getappdata(handles.leadfigure, 'M');

        M.patient.list=[M.patient.list; folders];
        M.patient.group=[M.patient.group; ones(length(folders),1)];

        setappdata(handles.leadfigure, 'M', M);
        ea_refresh_lg_connectome(handles);
        % save M
        M=getappdata(handles.leadfigure, 'M');
        save([handles.groupdir_choosebox.String, 'GroupConnectomeAnalysis.mat'], 'M', '-v7.3');
    end
end


% --- Outputs from this function are returned to the command line.
function varargout = lead_group_connectome_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in patientlist.
function patientlist_Callback(hObject, eventdata, handles)
% hObject    handle to patientlist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns patientlist contents as cell array
%        contents{get(hObject,'Value')} returns selected item from patientlist
M=getappdata(gcf,'M');

M.ui.listselect=get(handles.patientlist,'Value');
set(handles.grouplist,'Value',M.ui.listselect);

setappdata(gcf,'M',M);
ea_refresh_lg_connectome(handles);


% --- Executes during object creation, after setting all properties.
function patientlist_CreateFcn(hObject, eventdata, handles)
% hObject    handle to patientlist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in addptbutton.
function addptbutton_Callback(hObject, eventdata, handles)
% hObject    handle to addptbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if strcmp(handles.groupdir_choosebox.String, 'Choose Group Directory')
    ea_error('Please choose a group directory first to store the group analysis!', simpleStack = 1);
end

M=getappdata(handles.leadfigure,'M');

folders=ea_uigetdir(ea_startpath,'Select Patient folders..');
M.patient.list=[M.patient.list;folders'];
M.patient.group=[M.patient.group;ones(length(folders),1)];

setappdata(handles.leadfigure,'M',M);
ea_refresh_lg_connectome(handles);
% save M
M=getappdata(handles.leadfigure,'M');
save([handles.groupdir_choosebox.String,'GroupConnectomeAnalysis.mat'],'M','-v7.3');


% --- Executes on button press in removeptbutton.
function removeptbutton_Callback(hObject, eventdata, handles)
% hObject    handle to removeptbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
M=getappdata(gcf,'M');

deleteentry=get(handles.patientlist,'Value');

M.patient.list(deleteentry)=[];

M.patient.group(deleteentry)=[];

for cvar=1:length(M.clinical.vars)
    try
        M.clinical.vars{cvar}(deleteentry,:)=[];
    end
end

setappdata(gcf,'M',M);
ea_refresh_lg_connectome(handles);


% --- Executes on selection change in clinicallist.
function clinicallist_Callback(hObject, eventdata, handles)
% hObject    handle to clinicallist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns clinicallist contents as cell array
%        contents{get(hObject,'Value')} returns selected item from clinicallist
M=getappdata(gcf,'M');

M.ui.clinicallist=get(handles.clinicallist,'Value');
setappdata(gcf,'M',M);
ea_refresh_lg_connectome(handles);


% --- Executes during object creation, after setting all properties.
function clinicallist_CreateFcn(hObject, eventdata, handles)
% hObject    handle to clinicallist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in addvarbutton.
function addvarbutton_Callback(hObject, eventdata, handles)
% hObject    handle to addvarbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
M=getappdata(gcf,'M');
M.ui.clinicallist=length(M.clinical.labels)+1;
[numat,nuvar]=ea_get_clinical(M);
if ~isempty(numat) % user did not press cancel
    M.clinical.vars{end+1}=numat;
    M.clinical.labels{end+1}=nuvar;
end
set(handles.clinicallist,'Value',M.ui.clinicallist);
% store and refresh UI
setappdata(gcf,'M',M);

ea_refresh_lg_connectome(handles);


function [mat,matname]=ea_get_clinical(M)
try
    mat=M.clinical.vars{M.ui.clinicallist};
catch % new variable
    mat=[];
end
try
    matname=M.clinical.labels{M.ui.clinicallist};
catch
    matname='New variable';
end
[numat,nuname]=ea_edit_regressor(M);

if ~isempty(numat) % user did not press cancel
    mat=numat;
    matname=nuname;
end


% --- Executes on button press in removevarbutton.
function removevarbutton_Callback(hObject, eventdata, handles)
% hObject    handle to removevarbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
M=getappdata(gcf,'M');

% delete data
M.clinical.vars(get(handles.clinicallist,'Value'))=[];
M.clinical.labels(get(handles.clinicallist,'Value'))=[];

% store and refresh UI
setappdata(gcf,'M',M);
ea_refresh_lg_connectome(handles);


function [pathname] = ea_uigetdir(start_path, dialog_title)
% Pick a directory with the Java widgets instead of uigetdir

import javax.swing.JFileChooser;

if nargin == 0 || strcmp(start_path,'') % || start_path == 0 % Allow a null argument.
    start_path = pwd;
end

jchooser = javaObjectEDT('javax.swing.JFileChooser', start_path);

jchooser.setFileSelectionMode(JFileChooser.FILES_AND_DIRECTORIES);
if nargin > 1
    jchooser.setDialogTitle(dialog_title);
end

jchooser.setMultiSelectionEnabled(true);

status = jchooser.showOpenDialog([]);

if status == JFileChooser.APPROVE_OPTION
    jFile = jchooser.getSelectedFiles();
    pathname{size(jFile, 1)}=[];
    for i=1:size(jFile, 1)
        pathname{i} = char(jFile(i).getAbsolutePath);
    end

elseif status == JFileChooser.CANCEL_OPTION
    pathname = [];
else
    error('Error occured while picking file.');
end


% --- Executes on selection change in grouplist.
function grouplist_Callback(hObject, eventdata, handles)
% hObject    handle to grouplist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns grouplist contents as cell array
%        contents{get(hObject,'Value')} returns selected item from grouplist
M=getappdata(gcf,'M');

M.ui.listselect=get(handles.grouplist,'Value');

set(handles.patientlist,'Value',M.ui.listselect);

setappdata(gcf,'M',M);
ea_refresh_lg_connectome(handles);


% --- Executes during object creation, after setting all properties.
function grouplist_CreateFcn(hObject, eventdata, handles)
% hObject    handle to grouplist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in plusgroupbutton.
function plusgroupbutton_Callback(hObject, eventdata, handles)
% hObject    handle to plusgroupbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
M=getappdata(gcf,'M');
M.patient.group(get(handles.patientlist,'Value'))=M.patient.group(get(handles.patientlist,'Value'))+1;
setappdata(gcf,'M',M);
ea_refresh_lg_connectome(handles);


% --- Executes on button press in minusgroupbutton.
function minusgroupbutton_Callback(hObject, eventdata, handles)
% hObject    handle to minusgroupbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
M=getappdata(gcf,'M');
if M.patient.group(get(handles.patientlist,'Value'))>1
    M.patient.group(get(handles.patientlist,'Value'))=M.patient.group(get(handles.patientlist,'Value'))-1;
end
setappdata(gcf,'M',M);
ea_refresh_lg_connectome(handles);


% --- Executes on button press in reviewvarbutton.
function reviewvarbutton_Callback(hObject, eventdata, handles)
% hObject    handle to reviewvarbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
M=getappdata(gcf,'M');

% store variables
%M.clinical.vars{get(handles.clinicallist,'Value')}(isnan(M.clinical.vars{get(handles.clinicallist,'Value')}))=0;
[M.clinical.vars{get(handles.clinicallist,'Value')},M.clinical.labels{get(handles.clinicallist,'Value')}]=ea_get_clinical(M);


% store and refresh UI
setappdata(gcf,'M',M);

ea_refresh_lg_connectome(handles);


% --- Executes on button press in moveptdownbutton.
function moveptdownbutton_Callback(hObject, eventdata, handles)
% hObject    handle to moveptdownbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
M=getappdata(gcf,'M');
whichmoved=get(handles.patientlist,'Value');
if length(whichmoved)>1; return; end % more than one selected..
if whichmoved==1 % first entry anyways
    return
end
ix=1:length(M.patient.list);
ix(whichmoved)=ix(whichmoved)-1;
ix(whichmoved-1)=ix(whichmoved-1)+1;

M.patient.list=M.patient.list(ix);
M.patient.group=M.patient.group(ix);
setappdata(gcf,'M',M);
set(handles.patientlist,'Value',whichmoved-1);

ea_refresh_lg_connectome(handles);


% --- Executes on button press in moveptupbutton.
function moveptupbutton_Callback(hObject, eventdata, handles)
% hObject    handle to moveptupbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

M=getappdata(gcf,'M');
whichmoved=get(handles.patientlist,'Value');
if length(whichmoved)>1; return; end % more than one selected..
if whichmoved==length(M.patient.list) % last entry anyways
    return
end
ix=1:length(M.patient.list);
ix(whichmoved)=ix(whichmoved)+1;
ix(whichmoved+1)=ix(whichmoved+1)-1;

M.patient.list=M.patient.list(ix);
M.patient.group=M.patient.group(ix);
setappdata(gcf,'M',M);
set(handles.patientlist,'Value',whichmoved+1);

ea_refresh_lg_connectome(handles);


% --- Executes on selection change in fiberspopup.
function fiberspopup_Callback(hObject, eventdata, handles)
% hObject    handle to fiberspopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns fiberspopup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from fiberspopup
M=getappdata(gcf,'M');
M.ui.connectomename = eventdata.Source.String{eventdata.Source.Value};
setappdata(gcf,'M',M);
ea_refresh_lg_connectome(handles);


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


% --- Executes on selection change in labelpopup.
function labelpopup_Callback(hObject, eventdata, handles)
% hObject    handle to labelpopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns labelpopup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from labelpopup
M=getappdata(gcf,'M');
M.ui.labelpopup = eventdata.Source.String{eventdata.Source.Value};
setappdata(gcf,'M',M);
ea_refresh_lg_connectome(handles);


% --- Executes during object creation, after setting all properties.
function labelpopup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to labelpopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in groupdir_choosebox.
function groupdir_choosebox_Callback(hObject, eventdata, handles)
% hObject    handle to groupdir_choosebox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% groupdir=ea_uigetdir(ea_startpath,'Choose Group Directory');
groupdir = uigetdir;

if ~groupdir % user pressed cancel
    return
end

ea_busyaction('on',handles.leadfigure,'group');

groupdir=[groupdir,filesep];
M=ea_initializeM_connectome;

set(handles.groupdir_choosebox,'String',groupdir);

try % if file already exists, load it (and overwrite M).
    load([groupdir,'GroupConnectomeAnalysis.mat']);
catch % if not, store it saving M.
    save([groupdir,'GroupConnectomeAnalysis.mat'],'M','-v7.3');
end

M.root = groupdir;
setappdata(handles.leadfigure,'M',M);

ea_busyaction('off',handles.leadfigure,'group');

ea_refresh_lg_connectome(handles);


% --- Executes on button press in opensubgui.
function opensubgui_Callback(hObject, eventdata, handles)
% hObject    handle to opensubgui (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
M=getappdata(gcf,'M');
[selection]=ea_groupselectorwholelist(M.ui.listselect,M.patient.list);

lead_dbs('loadsubs',M.patient.list(selection));


% --- Executes on selection change in normregpopup.
function normregpopup_Callback(hObject, eventdata, handles)
% hObject    handle to normregpopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns normregpopup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from normregpopup
M=getappdata(gcf,'M');
M.ui.normregpopup=get(handles.normregpopup,'Value');
setappdata(gcf,'M',M);
ea_refresh_lg_connectome(handles);


% --- Executes during object creation, after setting all properties.
function normregpopup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to normregpopup (see GCBO)
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

ea_busyaction('on',gcf,'group');
if ~strcmp(handles.groupdir_choosebox.String,'Choose Group Directory') % group dir still not chosen
    disp('Saving data...');
    % save M
    ea_refresh_lg_connectome(handles);
    M=getappdata(hObject,'M');
    try
        save([handles.groupdir_choosebox.String,'GroupConnectomeAnalysis.mat'],'M','-v7.3');
    catch
        warning('Data could not be saved.');
        keyboard
    end
    disp('Done.');
    disp('Bye for now.');
end
ea_busyaction('off',gcf,'group');
delete(hObject);


% --- Executes on button press in lc_SPM.
function lc_SPM_Callback(hObject, eventdata, handles)
% hObject    handle to lc_SPM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
M=getappdata(handles.leadfigure,'M');

gecs = get(handles.lc_graphmetric,'String');
gecs = gecs{M.ui.lc.graphmetric};
parc = get(handles.labelpopup,'String');
parc = parc{get(handles.labelpopup,'Value')};
if M.ui.lc.smooth
    smoothsuffix='_smoothed';
else
    smoothsuffix='';
end


twosample=any(ismember(gecs,'>'));
if twosample
    cpair=ea_strsplit(gecs,'>');
end
for sub=1:length(M.patient.list)
    if M.ui.lc.normalization == 1 % no normalization
        normflag = '';
    else
        if M.ui.lc.normalization == 2 % z-Score
            normflag = 'z';
        elseif M.ui.lc.normalization == 3 % Albada 2008
            normflag = 'k';
        end
        if twosample
            ea_histnormalize([M.patient.list{sub},filesep,'connectomics',filesep,parc,filesep,'graph',filesep,cpair{1},'.nii'], normflag);
            ea_histnormalize([M.patient.list{sub},filesep,'connectomics',filesep,parc,filesep,'graph',filesep,cpair{2},'.nii'], normflag);
        else
            ea_histnormalize([M.patient.list{sub},filesep,'connectomics',filesep,parc,filesep,'graph',filesep,gecs,'.nii'], normflag);
        end
    end

    tn=ea_load_nii([M.patient.list{sub},filesep,'connectomics',filesep,parc,filesep,'graph',filesep,cpair{1},'.nii']);
   if any(tn.voxsize>2)
       if twosample
           ea_reslice_nii([M.patient.list{sub},filesep,'connectomics',filesep,parc,filesep,'graph',filesep,normflag,cpair{1},'.nii'],...
               [M.patient.list{sub},filesep,'connectomics',filesep,parc,filesep,'graph',filesep,'re',normflag,cpair{1},'.nii'],[2,2,2]);
           ea_reslice_nii([M.patient.list{sub},filesep,'connectomics',filesep,parc,filesep,'graph',filesep,normflag,cpair{2},'.nii'],...
               [M.patient.list{sub},filesep,'connectomics',filesep,parc,filesep,'graph',filesep,'re',normflag,cpair{2},'.nii'],[2,2,2]);
       else
           ea_reslice_nii([M.patient.list{sub},filesep,'connectomics',filesep,parc,filesep,'graph',filesep,normflag,gecs,'.nii'],...
               [M.patient.list{sub},filesep,'connectomics',filesep,parc,filesep,'graph',filesep,'re',normflag,gecs,'.nii'],[2,2,2]);
       end
       normflag=['re',normflag];
   end



    if M.ui.lc.smooth
        smoothflag = 's';
        if twosample
            ea_smooth([M.patient.list{sub},filesep,'connectomics',filesep,parc,filesep,'graph',filesep,normflag,cpair{1},'.nii']);
            ea_smooth([M.patient.list{sub},filesep,'connectomics',filesep,parc,filesep,'graph',filesep,normflag,cpair{2},'.nii']);
        else
            ea_smooth([M.patient.list{sub},filesep,'connectomics',filesep,parc,filesep,'graph',filesep,normflag,gecs,'.nii']);
        end
    else
        smoothflag = '';
    end
if twosample
    fis{sub,1} = [M.patient.list{sub},filesep,'connectomics',filesep,parc,filesep,'graph',filesep,smoothflag,normflag,cpair{1},'.nii'];
    fis{sub,2} = [M.patient.list{sub},filesep,'connectomics',filesep,parc,filesep,'graph',filesep,smoothflag,normflag,cpair{2},'.nii'];
else
    fis{sub} = [M.patient.list{sub},filesep,'connectomics',filesep,parc,filesep,'graph',filesep,smoothflag,normflag,gecs,'.nii'];
end

end

spmdir = [M.root,'connectomics',filesep,parc,filesep,'graph',filesep,normflag,gecs,smoothsuffix,filesep,'SPM'];
if exist(spmdir, 'dir')
    rmdir(spmdir,'s');
end
mkdir(spmdir);

%% model specification:
if twosample
    matlabbatch{1}.spm.stats.factorial_design.dir = {spmdir};
    for p=1:size(fis,1)
        matlabbatch{1}.spm.stats.factorial_design.des.pt.pair(p).scans = fis(p,:)';
    end
    matlabbatch{1}.spm.stats.factorial_design.des.pt.gmsca = 0;
    matlabbatch{1}.spm.stats.factorial_design.des.pt.ancova = 0;
    matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
    matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
    matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
    matlabbatch{1}.spm.stats.factorial_design.masking.im = 0;
    matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
    matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
    matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
    matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
else
    matlabbatch{1}.spm.stats.factorial_design.dir = {spmdir};
    matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = fis';
    matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
    matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
    matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
    matlabbatch{1}.spm.stats.factorial_design.masking.im = 0;
    matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
    matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
    matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
    matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
end
    spm_jobman('run',{matlabbatch});
    clear matlabbatch

%% model estimation:
matlabbatch{1}.spm.stats.fmri_est.spmmat = {[spmdir, filesep, 'SPM.mat']};
matlabbatch{1}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;
spm_jobman('run',{matlabbatch});
clear matlabbatch

%% contrast manager:
if twosample
    matlabbatch{1}.spm.stats.con.spmmat = {[spmdir, filesep, 'SPM.mat']};
    matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = gecs;
    matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = [1 -1];
    matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
    matlabbatch{1}.spm.stats.con.consess{2}.tcon.name = strrep(gecs,'>','<');
    matlabbatch{1}.spm.stats.con.consess{2}.tcon.weights = [-1 1];
    matlabbatch{1}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
    matlabbatch{1}.spm.stats.con.delete = 1;
else
    matlabbatch{1}.spm.stats.con.spmmat = {[spmdir, filesep, 'SPM.mat']};
    matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'main effect';
    matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = 1;
    matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
    matlabbatch{1}.spm.stats.con.delete = 1;
end
spm_jobman('run',{matlabbatch});
clear matlabbatch


% --- Executes on selection change in lc_normalization.
function lc_normalization_Callback(hObject, eventdata, handles)
% hObject    handle to lc_normalization (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns lc_normalization contents as cell array
%        contents{get(hObject,'Value')} returns selected item from lc_normalization
M=getappdata(gcf,'M');
M.ui.lc.normalization=get(handles.lc_normalization,'Value');
setappdata(gcf,'M',M);

% --- Executes during object creation, after setting all properties.
function lc_normalization_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lc_normalization (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in lc_graphmetric.
function lc_graphmetric_Callback(hObject, eventdata, handles)
% hObject    handle to lc_graphmetric (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns lc_graphmetric contents as cell array
%        contents{get(hObject,'Value')} returns selected item from lc_graphmetric
M=getappdata(gcf,'M');
M.ui.lc.graphmetric=get(handles.lc_graphmetric,'Value');
setappdata(gcf,'M',M);

% --- Executes during object creation, after setting all properties.
function lc_graphmetric_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lc_graphmetric (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in lc_smooth.
function lc_smooth_Callback(hObject, eventdata, handles)
% hObject    handle to lc_smooth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of lc_smooth
M=getappdata(gcf,'M');
M.ui.lc.smooth=get(handles.lc_smooth,'Value');
setappdata(gcf,'M',M);


function ea_histnormalize(fname, normflag)
nii=ea_load_nii(fname);
[pth,fn,ext]=fileparts(fname);

msk = isnan(nii.img) | (nii.img == 0);
vals = nii.img(~msk);

switch normflag
    case 'z' % zscore
        nii.fname=[pth,filesep,'z',fn,ext];
        vals=zscore(vals(:));
    case 'k' % albada
        nii.fname=[pth,filesep,'k',fn,ext];
        vals=ea_normal(vals(:));
end

nii.img(~msk)=vals;
spm_write_vol(nii,nii.img);


function ea_smooth(fname)
[pth,fn,ext]=fileparts(fname);
spm_smooth(fname,[pth,filesep,'s',fn,ext],[8 8 8]);


% --- Executes on button press in calcgroupconnectome.
function calcgroupconnectome_Callback(hObject, eventdata, handles)
% hObject    handle to calcgroupconnectome (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

options.prefs=ea_prefs('tmp');
M=getappdata(gcf,'M');

disp('Concatenating connectome:');
howmanyfibs=inputdlg('How many fibers to sample from each subject?','Sample fibers',1,{'20000'});
howmanyfibs=str2double(howmanyfibs);

ftrFiles = cell(length(M.patient.list), 1);
for sub=1:length(M.patient.list)
    ftrFiles{sub} = [M.patient.list{sub},filesep,'connectomes',filesep,'dMRI',filesep,options.prefs.FTR_normalized];
end

if ~exist([M.root,'connectomes',filesep,'dMRI'], 'dir')
    mkdir([M.root,'connectomes',filesep,'dMRI'])
end

entries=get(handles.gcfilter,'String');
entry=entries{get(handles.gcfilter,'Value')};
switch entry
    case 'No filtering'
        filtermask='';
    case 'Filter by white matter mask'
        filtermask=ea_niigz([ea_space,'c2mask']);
    case 'Filter by brain mask'
        filtermask=ea_niigz([ea_space,'brainmask']);
    case 'Filter by custom nii...'
        filtermask=getappdata(handles.gcfilter,'filtermask');
end

ea_ftr_aggregate(ftrFiles, ...
    [M.root,'connectomes',filesep,'dMRI',filesep,options.prefs.FTR_normalized], ...
    howmanyfibs, 'number', filtermask);


function lc_contrast_Callback(hObject, eventdata, handles)
% hObject    handle to lc_contrast (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of lc_contrast as text
%        str2double(get(hObject,'String')) returns contents of lc_contrast as a double


% --- Executes during object creation, after setting all properties.
function lc_contrast_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lc_contrast (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in lc_stattest.
function lc_stattest_Callback(hObject, eventdata, handles)
% hObject    handle to lc_stattest (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns lc_stattest contents as cell array
%        contents{get(hObject,'Value')} returns selected item from lc_stattest


% --- Executes during object creation, after setting all properties.
function lc_stattest_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lc_stattest (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in lc_metric.
function lc_metric_Callback(hObject, eventdata, handles)
% hObject    handle to lc_metric (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns lc_metric contents as cell array
%        contents{get(hObject,'Value')} returns selected item from lc_metric


% --- Executes during object creation, after setting all properties.
function lc_metric_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lc_metric (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function lc_threshold_Callback(hObject, eventdata, handles)
% hObject    handle to lc_threshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of lc_threshold as text
%        str2double(get(hObject,'String')) returns contents of lc_threshold as a double


% --- Executes during object creation, after setting all properties.
function lc_threshold_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lc_threshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in lc_nbs.
function lc_nbs_Callback(hObject, eventdata, handles)
% hObject    handle to lc_nbs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

clearvars -global nbs
global nbs

earoot=getappdata(handles.leadfigure,'earoot');
UI.method.ui = 'Run NBS';
UI.test.ui = get(handles.lc_stattest,'String');
UI.test.ui=UI.test.ui{get(handles.lc_stattest,'Value')};
UI.thresh.ui = get(handles.lc_threshold,'String');
UI.contrast.ui = get(handles.lc_contrast,'String');
ea_preparenbs(handles)
root=handles.groupdir_choosebox.String;
UI.design.ui = [root,'NBSdesignMatrix.mat'];
UI.matrices.ui = [root,'NBSdataMatrix.mat'];
UI.node_coor.ui = '';
UI.node_label.ui = '';

% get advanced options:
try
    lc=load([earoot,'connectomics',filesep,'lc_options.mat']);
catch
    lc=ea_initlcopts([]);
end
if ~isfield(lc,'nbs') % compatibility with older stored userdata (<v1.4.9)
    lc=ea_initlcopts([],lc); % will merely add the nbs stuff
end
save([earoot,'connectomics',filesep,'lc_options.mat'],'-struct','lc');
switch lc.nbs.adv.compsize
    case 1
        UI.size.ui = 'Extent';
    case 2
        UI.size.ui = 'Intensity';
end

UI.perms.ui = num2str(lc.nbs.adv.perm);
UI.alpha.ui = num2str(lc.nbs.adv.alpha);
UI.exchange.ui = lc.nbs.adv.exch;

ea_NBSrun(UI,[]);

load([root,'NBSdataMatrix.mat']);
load([root,'NBSdesignMatrix.mat']);

switch UI.test.ui
    case 't-test'
        c=eval(UI.contrast.ui);
        if length(c)>2
            ea_error('Only two-sample t-tests are fully supported at present.');
        end

        [~,~,~,tstat]=ttest2(permute(allX(:,:,logical(mX(:,find(c==1)))),[3,1,2]),...
            permute(allX(:,:,logical(mX(:,find(c==-1)))),[3,1,2]));

        T=squeeze(tstat.tstat);
        uT=T;
        clear tstat
        clear pmask
        pmask=zeros([nbs.NBS.n,size(T)]);
        for network=1:nbs.NBS.n
            X=full(nbs.NBS.con_mat{network});
            X=X+X';
            pmask(network,:,:)=X;
            save([root,'sig_',num2str(network)],'X');
        end

        if network
            pmask=squeeze(sum(pmask,1));
        end
        T(~pmask)=nan;
    otherwise

        for network=1:nbs.NBS.n
            X=full(nbs.NBS.con_mat{network});
            X=X+X';
            save([root,'sig_',num2str(network)],'X');
        end
        ea_error('Only t-test fully supported at present. Please process results manually for other tests.');
end

thisparc=get(handles.labelpopup,'String');
thisparc=thisparc{get(handles.labelpopup,'Value')};
thismetr=get(handles.lc_metric,'String');
thismetr=thismetr{get(handles.lc_metric,'Value')};

if ismember('&',thismetr) % clean from & for variable
    [~,ix]=ismember('&',thismetr);
    thismetr(ix)='';
end
eval([thismetr,'=T;']);
expfn=[root,thisparc,'_',thismetr,'T_p<',UI.alpha.ui,'.mat'];
save(expfn,thismetr);
eval([thismetr,'=uT;']);
expfn=[root,thisparc,'_',thismetr,'T_unthresholded.mat'];
save(expfn,thismetr);

disp('** NBS done.');
if network
    disp('NBS found at least one significant network.');
    disp(['It has been stored in: ',expfn,'.']);

    disp(['To display the network(s), please load this file in the 3D viewer''s "Connectivity Visualization" under the "Matrix Level" panel.']);
else
    disp('NBS found no significant networks.');
end


function ea_preparenbs(handles)

% prepare designmatrix:

gstr=(get(handles.grouplist,'String'));

for pt=1:length(gstr)
    gv(pt)=str2double(gstr(pt));
end
mX=zeros(length(gv),max(gv));
for g=1:max(gv)
    mX(:,g)=gv==g;
end
save([handles.groupdir_choosebox.String,'NBSdesignMatrix'],'mX');

% prepare data matrix:

M=getappdata(handles.leadfigure,'M');
thisparc=get(handles.labelpopup,'String');
thisparc=thisparc{get(handles.labelpopup,'Value')};

thismetr=get(handles.lc_metric,'String');
thismetr=thismetr{get(handles.lc_metric,'Value')};

if ismember('&',thismetr)
    % e.g. ON vs OFF metric - need to load two matrices per subject!
    [~,plusix]=ismember('&',thismetr);
    [scix]=strfind(thismetr,'_');
    suffx1=thismetr(scix(1)+1:plusix-1);
    suffx2=thismetr(plusix+1:scix(2)-1);
    base=thismetr(1:scix(1)-1);
    cap=thismetr(scix(2)+1:end);
    clear thismetr
    metrix{1}=[base,'_',suffx1,'_',cap];
    metrix{2}=[base,'_',suffx2,'_',cap];

    % inflate design matrix
    if ~all(mX==1)
        ea_error('Multi subject measurements only supported for a single cohort. All group values must be 1.');
    end

    mX=repmat(eye(2),length(mX),1);
    save([handles.groupdir_choosebox.String,'NBSdesignMatrix'],'mX');
else
 metrix{1}=thismetr;
end

Xcnt=1;
for pt=1:length(M.patient.list)
    for metr=1:length(metrix)
        X=load([M.patient.list{pt},filesep,'connectomics',filesep,thisparc,filesep,metrix{metr},'.mat']);
        fn=fieldnames(X);
        if ~exist('allX','var')
            allX=nan([size(X.(fn{1})),length(M.patient.list)]);
        end

        X=X.(fn{1});
        switch get(handles.normregpopup,'Value')
            case 2
                X(:)=ea_nanzscore(X(:));
            case 3
                X(:)=ea_normal(X(:));
        end

        allX(:,:,Xcnt)=X;
        Xcnt=Xcnt+1;
    end
end
save([handles.groupdir_choosebox.String,'NBSdataMatrix'],'allX','-v7.3');


% --- Executes on button press in lc_nbsadvanced.
function lc_nbsadvanced_Callback(hObject, eventdata, handles)
% hObject    handle to lc_nbsadvanced (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

ea_nbs_advanced;


% --- Executes on selection change in gcfilter.
function gcfilter_Callback(hObject, eventdata, handles)
% hObject    handle to gcfilter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns gcfilter contents as cell array
%        contents{get(hObject,'Value')} returns selected item from gcfilter

entries=get(hObject,'String');
entry=entries{get(hObject,'Value')};
if strcmp(entry,'Filter by custom nii...')
    [file,path]=uigetfile({'*.nii','.nii.gz'},'Select nifti file for filtering...');
    setappdata(hObject,'filtermask',fullfile(path,file));
end

% --- Executes during object creation, after setting all properties.
function gcfilter_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gcfilter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
