function varargout = lead(varargin)
% LEAD M-file for lead.fig
%      LEAD, by itself, creates a new LEAD or raises the existing
%      singleton*.
%
%      H = LEAD returns the handle to a new LEAD or the handle to
%      the existing singleton*.
%
%      LEAD('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in LEAD.M with the given input arguments.
%
%      LEAD('Property','Value',...) creates a new LEAD or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before lead_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to lead_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help lead

% Last Modified by GUIDE v2.5 29-Feb-2016 16:44:57

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @lead_OpeningFcn, ...
                   'gui_OutputFcn',  @lead_OutputFcn, ...
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


% --- Executes just before lead is made visible.
function lead_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to lead (see VARARGIN)

% Choose default command line output for lead
handles.output = hObject;


% Update handles structure
guidata(hObject, handles);
addpath(genpath(fileparts(which(mfilename))));
ea_dispbn;


% load atlassets
earoot=[fileparts(which('lead')),filesep];
as=dir([earoot,'atlases',filesep]);
asc=cell(0);
cnt=1;
for i=1:length(as)
    if as(i).isdir
    asc{cnt}=as(i).name;
    cnt=cnt+1;
    end
end
options.prefs=ea_prefs('');
if options.prefs.env.dev
asc{end+1}='Segment patient anatomy';
end
excludes={'.','..'};
asc(ismember(asc,excludes))=[];

asc{end+1}='Use none';

set(handles.atlassetpopup,'String',asc);


set(handles.normalize_checkbox,'Value',0);

set(hObject,'Color',[1 1 1]);
set(handles.versiontxt,'String',['v',ea_getvsn('local')]);


% set DICOM input and output name strings:
try
load([fileparts(which('lead')),filesep,'ea_prefs']);
set(handles.setdicomout,'String',lp.dicom.outfolder);
set(handles.setdicomin,'String',lp.dicom.infolder);
end


% check if group connectome files are present
if ~exist([earoot,'fibers'],'file') || ...
        ~exist([earoot,'fibers',filesep,'Groupconnectome (Horn 2013) full.mat'],'file') || ...
        ~exist([earoot,'fibers',filesep,'Groupconnectome (Horn 2013) thinned out x 2.mat'],'file') || ...
        ~exist([earoot,'fibers',filesep,'Groupconnectome (Horn 2013) thinned out x 5.mat'],'file') || ...
        ~exist([earoot,'fibers',filesep,'Groupconnectome (Horn 2013) thinned out x 5.mat'],'file') || ...
        ~exist([earoot,'fibers',filesep,'Groupconnectome (Horn 2013) thinned out x 10.mat'],'file') || ...
        ~exist([earoot,'fibers',filesep,'Groupconnectome (Horn 2013) thinned out x 50.mat'],'file') || ...
        ~exist([earoot,'fibers',filesep,'Groupconnectome (Horn 2013) thinned out x 100.mat'],'file') || ...
        ~exist([earoot,'fibers',filesep,'Groupconnectome (Horn 2013) thinned out x 500.mat'],'file')
    set(handles.dlgroupc,'Visible','on');
    set(handles.dlgroupc,'BackgroundColor',[0.2,0.8,0.2]);
else
    set(handles.dlgroupc,'Visible','off');
end


%im = imread('bg_gui.png');
%image(im);
%axis off;
%axis fill

%set(0,'gca',handles.logoaxes);
im=imread('ea_logo.png');
image(im);
axis off;
axis equal;

try
warning('off');
set(handles.dicompanel,'BackgroundColor','none');
set(handles.gopanel,'BackgroundColor','none');
set(handles.psapanel,'BackgroundColor','none');
set(handles.normpanel,'BackgroundColor','none');
set(handles.reconpanel,'BackgroundColor','none');
set(handles.reviewpanel,'BackgroundColor','none');
set(handles.vizpanel,'BackgroundColor','none');
set(handles.coregctpanel,'BackgroundColor','none');
warning('on');
end

% get electrode model specs and place in popup
set(handles.electrode_model_popup,'String',ea_resolve_elspec);

set(gcf,'name','Welcome to LEAD-DBS');

% add normalization methods to menu
cnt=1;
ndir=dir([earoot,'ea_normalize_*.m']);

for nd=length(ndir):-1:1
    [~,methodf]=fileparts(ndir(nd).name);
    try
        [thisndc,spmvers]=eval([methodf,'(','''prompt''',')']);
        if ismember(spm('ver'),spmvers)
        ndc{cnt}=thisndc;
        normmethod{cnt}=methodf;
        if strcmp(ndc{cnt},eval([options.prefs.normalize.default,'(','''prompt''',')']))
         defentry=cnt;   
        end
        cnt=cnt+1;
        end
    end
end



try
setappdata(gcf,'normmethod',normmethod);
set(handles.normmethod,'String',ndc);
catch
    if isempty(which('spm'))
    warning('It seems that SPM is not installed.');
    end
end
try % set selection of normmethod to default entry (specified in ea_prefs).
    if defentry<=length(get(handles.normmethod,'String'))
        set(handles.normmethod,'Value',defentry);
    end
end
clear defentry

% add coreg methods to menu
cnt=1;
ndir=dir([earoot,'ea_coregctmri_*.m']);
for nd=length(ndir):-1:1
    [~,methodf]=fileparts(ndir(nd).name);
    try
        [thisndc,spmvers]=eval([methodf,'(','''prompt''',')']);
        if ismember(spm('ver'),spmvers)
        cdc{cnt}=thisndc;
        coregctmethod{cnt}=methodf;
        if strcmp(cdc{cnt},eval([options.prefs.ctcoreg.default,'(','''prompt''',')']))
         defentry=cnt;   
        end
        cnt=cnt+1;
        end
    end
end
try
setappdata(gcf,'coregctmethod',coregctmethod);
set(handles.coregctmethod,'String',cdc);
catch
    if isempty(which('spm'))
    warning('It seems that SPM is not installed.');
    end
end
try % set selection of ctcoregmethod to default entry (specified in ea_prefs).
    if defentry<=length(get(handles.coregctmethod,'String'))
        set(handles.coregctmethod,'Value',defentry);
    end
end



if nargin
    
    if ~isempty(varargin)
        switch varargin{1}
            case 'loadsubs'
                
                ea_load_pts(handles,varargin{2});
                
        end
    end
    
end


%% add tools menu
menuprobe=getappdata(handles.leadfigure,'menuprobe');
if isempty(menuprobe)
f = uimenu('Label','Tools');
    uimenu(f,'Label','Convert ACPC/MNI coordinates','Callback',{@ea_acpcquery,handles.leadfigure});
setappdata(handles.leadfigure,'menuprobe',1);
end

ea_firstrun(handles);
getui(handles);




% UIWAIT makes lead wait for user response (see UIRESUME)
% uiwait(handles.leadfigure);


% --- Outputs from this function are returned to the command line.
function varargout = lead_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in run_button.
function run_button_Callback(hObject, eventdata, handles)
% hObject    handle to run_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% run trajectory reconstruction.
leadfig=handles.leadfigure;
ea_busyaction('on',leadfig,'lead');

options=handles2options(handles);

try
    options.lc=load([fileparts(which('lead')),filesep,'connectomics',filesep,'lc_options.mat']);
catch
    options.lc=[];
end

if options.d3.autoserver && options.d3.write
    choice = questdlg('Are you sure you want to export results to the server automatically?', ...
        'Auto Server-Export', ...
        'Cancel','Yes','Cancel');
    % Handle response
    switch choice
        case 'Cancel'
            return
    end
end

clc
uipatdirs=getappdata(gcf,'uipatdir');

if isempty(uipatdirs)
    uipatdirs={'No Patient Selected'};
end

prefs=ea_prefs('');
if length(uipatdirs)>1 && ~isempty(which('parpool')) && prefs.pp.do % do parallel processing if available and set in ea_prefs.
try delete(gcp); end
    pp=parpool(prefs.pp.profile,prefs.pp.csize);

    for pat=1:length(uipatdirs)
        % set patient specific options
        opts{pat}=options;
        opts{pat}.root=[fileparts(uipatdirs{pat}),filesep];
        [~,thispatdir]=fileparts(uipatdirs{pat});
        opts{pat}.patientname=thispatdir;
    end

    parfor pat=1:length(uipatdirs)

        % run main function
        try
            ea_autocoord(opts{pat});
        catch
            warning([opts{pat}.patientname,' failed. Please run this patient again and adjust parameters. Moving on to next patient.' ]);
        end

    end
    delete(pp);

else

    for pat=1:length(uipatdirs)
        % set patient specific options
        options.root=[fileparts(uipatdirs{pat}),filesep];
        [root,thispatdir]=fileparts(uipatdirs{pat});
        options.patientname=thispatdir;
        % run main function

        if length(uipatdirs)>1 % multi mode. Dont stop at errors.
            try
                ea_autocoord(options);
            catch
                warning([options.patientname,' failed. Please run this patient again and adjust parameters. Moving on to next patient.' ]);
            end
        else
            ea_autocoord(options);
        end
    end
end
ea_busyaction('off',leadfig,'lead');




function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double
storeui(handles);


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
storeui(handles);



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
storeui(handles);


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
storeui(handles);


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
storeui(handles);


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
storeui(handles);


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
storeui(handles);


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
storeui(handles);



function z_contrast_txt_Callback(hObject, eventdata, handles)
% hObject    handle to z_contrast_txt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of z_contrast_txt as text
%        str2double(get(hObject,'String')) returns contents of z_contrast_txt as a double
storeui(handles);


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
storeui(handles);


% --- Executes on button press in render_checkbox.
function render_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to render_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of render_checkbox
storeui(handles);


% --- Executes on button press in showatlases_checkbox.
function showatlases_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to showatlases_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of showatlases_checkbox
storeui(handles);


% --- Executes on button press in showfibers_checkbox.
function showfibers_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to showfibers_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of showfibers_checkbox
storeui(handles);


% --- Executes on button press in showconnectivities_checkbox.
function showconnectivities_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to showconnectivities_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of showconnectivities_checkbox
storeui(handles);


% --- Executes on button press in writeout2d_checkbox.
function writeout2d_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to writeout2d_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of writeout2d_checkbox
storeui(handles);


% --- Executes on button press in showatlases2d_checkbox.
function showatlases2d_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to showatlases2d_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of showatlases2d_checkbox
storeui(handles);


% --- Executes on button press in manualheight_checkbox.
function manualheight_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to manualheight_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of manualheight_checkbox
storeui(handles);


% --- Executes on button press in patdir_choosebox.
function patdir_choosebox_Callback(hObject, eventdata, handles)
% hObject    handle to patdir_choosebox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
p='/';
try
load([fileparts(which('lead')),filesep,'ea_prefs']);
p=lp.dicom.outfolder;
end

uipatdir=ea_uigetdir(p,'Please choose patient folder(s)...');

if isempty(uipatdir)
    return
end

ea_load_pts(handles,uipatdir);



function ea_load_pts(handles,uipatdir)


if length(uipatdir)>1
    set(handles.patdir_choosebox,'String',['Multiple (',num2str(length(uipatdir)),')']);
    set(handles.patdir_choosebox,'TooltipString',ea_strjoin(uipatdir,', '));
else
    set(handles.patdir_choosebox,'String',uipatdir{1});
    set(handles.patdir_choosebox,'TooltipString',uipatdir{1});
end

% store patient directories in figure


setappdata(gcf,'uipatdir',uipatdir);
ea_switchctmr(handles);

getui(handles); % update ui from patient
storeui(handles); % save in pt folder

% --- Executes on button press in left_checkbox.
function left_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to left_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of left_checkbox
storeui(handles);


% --- Executes on button press in right_checkbox.
function right_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to right_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of right_checkbox
storeui(handles);


% --- Executes on button press in normalize_checkbox.
function normalize_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to normalize_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of normalize_checkbox
storeui(handles);



function refinesteps_txt_Callback(~, eventdata, handles)
% hObject    handle to refinesteps_txt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of refinesteps_txt as text
%        str2double(get(hObject,'String')) returns contents of refinesteps_txt as a double
storeui(handles);


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
storeui(handles);


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

if get(hObject,'Value')


   set(handles.maskwindow_txt,'Enable','on');
set(handles.targetpopup,'Enable','on');



else

   set(handles.targetpopup,'Enable','off');
   set(handles.maskwindow_txt,'Enable','off');
end
storeui(handles);



% --- Executes on selection change in MRCT.
function MRCT_Callback(hObject, eventdata, handles)
% hObject    handle to MRCT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns MRCT contents as cell array
%        contents{get(hObject,'Value')} returns selected item from MRCT

ea_switchctmr(handles,get(hObject,'Value'));


storeui(handles);


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
storeui(handles);


% --- Executes on button press in cg25check.
function cg25check_Callback(hObject, eventdata, handles)
% hObject    handle to cg25check (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cg25check
storeui(handles);


% --- Executes on selection change in normmethod.
function normmethod_Callback(hObject, eventdata, handles)
% hObject    handle to normmethod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns normmethod contents as cell array
%        contents{get(hObject,'Value')} returns selected item from normmethod
storeui(handles);


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
storeui(handles);


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
storeui(handles);


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
storeui(handles);


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
storeui(handles);


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
storeui(handles);










% --- Executes on button press in cmappushbutton.
function cmappushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to cmappushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
colormapeditor;
storeui(handles);



% --- Executes on button press in canatlcheck.
function canatlcheck_Callback(hObject, eventdata, handles)
% hObject    handle to canatlcheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of canatlcheck
storeui(handles);


% --- Executes on button press in patatlcheck.
function patatlcheck_Callback(hObject, eventdata, handles)
% hObject    handle to patatlcheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of patatlcheck
storeui(handles);


% --- Executes on button press in normptatlascheck.
function normptatlascheck_Callback(hObject, eventdata, handles)
% hObject    handle to normptatlascheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of normptatlascheck
storeui(handles);


% --- Executes on button press in dicomcheck.
function dicomcheck_Callback(hObject, eventdata, handles)
% hObject    handle to dicomcheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of dicomcheck
storeui(handles);


% --- Executes on button press in setdicomin.
function setdicomin_Callback(hObject, eventdata, handles)
% hObject    handle to setdicomin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
storeui(handles);

p='';
try
load([fileparts(which('lead')),filesep,'ea_prefs']);
p=lp.dicom.infolder;
end

dicindir=uigetdir(p);

if ~dicindir
    return
end

set(handles.setdicomin,'String',dicindir);

try
load([fileparts(which('lead')),filesep,'ea_prefs']);
end
lp.dicom.infolder=[dicindir,filesep];
save([fileparts(which('lead')),filesep,'ea_prefs'],'lp');


% --- Executes on button press in setdicomout.
function setdicomout_Callback(hObject, eventdata, handles)
% hObject    handle to setdicomout (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
storeui(handles);

p='';
try
load([fileparts(which('lead')),filesep,'ea_prefs']);
p=lp.dicom.outfolder;
end

dicoutdir=uigetdir(p);

if ~dicoutdir
    return
end
set(handles.setdicomout,'String',dicoutdir);
try
load([fileparts(which('lead')),filesep,'ea_prefs']);
end
lp.dicom.outfolder=[dicoutdir,filesep];
save([fileparts(which('lead')),filesep,'ea_prefs'],'lp');



function [pathname] = ea_uigetdir(start_path, dialog_title)
% Pick a directory with the Java widgets instead of uigetdir

import javax.swing.JFileChooser;

if nargin == 0 || strcmp(start_path,'') % Allow a null argument.
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


% --- Executes on button press in genptatlascheck.
function genptatlascheck_Callback(hObject, eventdata, handles)
% hObject    handle to genptatlascheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of genptatlascheck
storeui(handles);


function storeui(handles)

try
chooseboxname=get(handles.patdir_choosebox,'String');
catch
    return
end
% determine if patientfolder is set
switch chooseboxname
    case 'Choose Patient Directory'
        outdir=[fileparts(which('lead')),filesep];
    otherwise
        if strcmp(chooseboxname(1:8),'Multiple')
                    outdir=[fileparts(which('lead')),filesep];

        else
        outdir=[get(handles.patdir_choosebox,'String'),filesep];
        end
end

updatestatus(handles);

options=handles2options(handles);
try save([outdir,'ea_ui'],'-struct','options'); end



function updatestatus(handles)
try
uipatdir=getappdata(gcf,'uipatdir');

set(handles.statusone,'String','One or more MR-/CT-volumes missing.');
modality=ea_checkctmrpresent(handles);

% check if MRCT popup is set correctly

if any(modality)

   if ~modality(get(handles.MRCT,'Value'))
       set(handles.MRCT,'ForegroundColor',[0.8,0.5,0.5]);
   else
       set(handles.MRCT,'ForegroundColor',[0.5,0.8,0.5]);
   end

end


% check for reconstructions
if exist([uipatdir{1},filesep,'ea_coords.fcsv'],'file') && exist([uipatdir{1},filesep,'ea_reconstruction.mat'],'file')
    set(handles.statustwo,'String','Fiducials and Trajectory information present in folder. Will be overwritten if "Reconstruct" is set.');
elseif exist([uipatdir{1},filesep,'ea_coords.fcsv'],'file') && ~exist([uipatdir{1},filesep,'ea_reconstruction.mat'],'file')
    set(handles.statustwo,'String','Fiducials information present in folder. Will be overwritten if "Reconstruct" is set.');
elseif ~exist([uipatdir{1},filesep,'ea_coords.fcsv'],'file') && exist([uipatdir{1},filesep,'ea_reconstruction.mat'],'file')
    set(handles.statustwo,'String','Trajectory information present in folder. Will be overwritten if "Reconstruct" is set.');
elseif ~exist([uipatdir{1},filesep,'ea_coords.fcsv'],'file') && ~exist([uipatdir{1},filesep,'ea_reconstruction.mat'],'file')
    set(handles.statustwo,'String','No reconstruction available in folder. Set "Reconstruct" to start.');
end
end




function ea_switchctmr(varargin)
% varargin: handles, switchto (1=MR, 2=CT; optional).
% being called when new patient is loaded.
handles=varargin{1};
switchto=0;

if nargin==1 % autodetect
    % check if MRCT popup is set correctly
modality=ea_checkctmrpresent(handles);
switchto=find(modality);
if any(modality)
   if ~modality(get(handles.MRCT,'Value'))
       set(handles.MRCT,'Value',switchto);
   else
   end
end
else
    % switch MRCT popup to specified handle.
    switchto=varargin{2};
end

if ~(sum(switchto>0)>1) && ~isempty(switchto) % e.g. MR and CT present
    switch switchto
        case 1 % MR
            set(handles.coregct_checkbox,'Enable','off');
            set(handles.coregct_checkbox,'Value',0);
            set(handles.coregctmethod,'Enable','off');
            set(handles.coregctcheck,'Enable','off');
            set(handles.coregctcheck,'Value',0);
            set(handles.coregthreshs,'Enable','off');
            set(handles.coregmrpopup,'Enable','on');
        case 2 % CT
            set(handles.coregct_checkbox,'Enable','on');
            set(handles.coregctmethod,'Enable','on');
            set(handles.coregctcheck,'Enable','on');
            set(handles.coregthreshs,'Enable','on');
            set(handles.coregmrpopup,'Enable','off');
    end
end
updatestatus(handles);

function getui(handles)



% determine if patientfolder is set
switch get(handles.patdir_choosebox,'String')
    case {'Choose Patient Directory','Multiple'}
        outdir=[fileparts(which('lead')),filesep];
    otherwise
        outdir=get(handles.patdir_choosebox,'String');
end
try

options=load([outdir,'ea_ui']);
options2handles(options,handles); % update UI
end



function options=handles2options(handles)

%% some manual options that can be set:


options.endtolerance=10; % how many slices to use with zero signal until end of electrode estimate.
options.sprungwert=4; % how far electrode centroid may be (in xy axis) from last to current slice.
options.refinesteps=0; % how often to re-iterate to reconstruct trajectory. More than 2 should usually not be beneficial. Use 0 to use the direct measurement.
options.tra_stdfactor=0.9; % Default: 0.9 - the lower this factor, the lower the threshold (more included pixels in tra process).
options.cor_stdfactor=1.0; % Default: 1.0 - the higher this factor, the lower the threshold (more included pixels in cor process).




%% set options

%uipatdir=get(handles.patdir_choosebox,'String');

options.earoot=[fileparts(which('lead')),filesep];
options.dicomimp=get(handles.dicomcheck,'Value');

options.normalize.do=(get(handles.normalize_checkbox,'Value') == get(handles.normalize_checkbox,'Max'));
options.normalize.method=getappdata(gcf,'normmethod');
options.normalize.method=options.normalize.method{get(handles.normmethod,'Value')};
options.normalize.methodn=get(handles.normmethod,'Value');

options.normalize.check=(get(handles.normcheck,'Value') == get(handles.normcheck,'Max'));

% coreg CT
options.coregct.do=(get(handles.coregct_checkbox,'Value') == get(handles.coregct_checkbox,'Max'));
options.coregct.method=getappdata(gcf,'coregctmethod');
options.coregct.method=options.coregct.method{get(handles.coregctmethod,'Value')};
options.coregct.methodn=get(handles.coregctmethod,'Value');
options.coregct.coregthreshs= eval( [ '[', get(handles.coregthreshs,'String'), ']' ] );

options.coregctcheck=get(handles.coregctcheck,'Value');


options.coregmr.method=get(handles.coregmrpopup,'Value');

% set modality (MR/CT) in options
options.modality = get(handles.MRCT,'Value');




options.verbose=3; % 4: Show figures but close them 3: Show all but close all figs except resultfig 2: Show all and leave figs open, 1: Show displays only, 0: Show no feedback.

%sidelog=[get(handles.right_checkbox,'Value') == get(handles.right_checkbox,'Max'),get(handles.left_checkbox,'Value') == get(handles.left_checkbox,'Max')];
%sidepos=[1,2];

%options.sides=sidepos(logical(sidelog)); %side=1 -> left electrode, side=2 -> right electrode. both: [1:2]
options.sides=1:2;

options.doreconstruction=(get(handles.doreconstruction_checkbox,'Value') == get(handles.doreconstruction_checkbox,'Max'));
if strcmp(get(handles.maskwindow_txt,'String'),'auto')
options.maskwindow=10; % initialize at 10
options.automask=1; % set automask flag
else
options.maskwindow=str2num(get(handles.maskwindow_txt,'String')); % size of the window that follows the trajectory
options.automask=0; % unset automask flag
end
options.autoimprove=0; % if true, templates will be modified.
options.axiscontrast=8; % if 8: use tra only but smooth it before. % if 9: use mean of cor and tra but smooth it. % if 10: use raw tra only.
options.zresolution=10; % voxels are being parcellated into this amount of portions.

options.atl.genpt=get(handles.vizspacepopup,'Value')==2; % generate patient specific atlases
options.atl.normalize=0; % normalize patient specific atlasset. This is not done anymore for now.
options.atl.can=get(handles.vizspacepopup,'Value')==1; % display canonical atlases
options.atl.pt=0; % display patient specific atlases. This is not done anymore for now.
options.atl.ptnative=get(handles.vizspacepopup,'Value')==2; % show results in native space.
if options.atl.ptnative
    options.native=1;
else
    options.native=0;
end

options.d2.write=(get(handles.writeout2d_checkbox,'Value') == get(handles.writeout2d_checkbox,'Max'));
options.d2.atlasopacity=0.15;


options.manualheightcorrection=(get(handles.manualheight_checkbox,'Value') == get(handles.manualheight_checkbox,'Max'));
options.d3.write=(get(handles.render_checkbox,'Value') == get(handles.render_checkbox,'Max'));
options.d3.prolong_electrode=2;
options.d3.verbose='on';
options.d3.elrendering=1;
options.d3.hlactivecontacts=0;
options.d3.showactivecontacts=1;
options.d3.showpassivecontacts=1;
options.d3.showisovolume=0;
options.d3.isovscloud=0;
options.d3.autoserver=get(handles.exportservercheck,'Value');

options.numcontacts=4;
options.entrypoint=get(handles.targetpopup,'String');
options.entrypoint=options.entrypoint{get(handles.targetpopup,'Value')};
options.entrypointn=get(handles.targetpopup,'Value');

options.writeoutpm=1;

options.elmodeln = get(handles.electrode_model_popup,'Value');
string_list = get(handles.electrode_model_popup,'String');
options.elmodel=string_list{options.elmodeln};
options.atlasset=get(handles.atlassetpopup,'String'); %{get(handles.atlassetpopup,'Value')}
options.atlasset=options.atlasset{get(handles.atlassetpopup,'Value')};
options.atlassetn=get(handles.atlassetpopup,'Value');

if strcmp(options.atlasset,'Use none');
    options.d3.writeatlases=0;
    options.d2.writeatlases=1;
else
    options.d3.writeatlases=1;
    options.d2.writeatlases=1;
end


options.expstatvat.do=0;

options.fiberthresh=10;
options.writeoutstats=1;


options.colormap=colormap;

options.dolc=get(handles.include_lead_connectome_subroutine,'Value');





function options2handles(options,handles)


%% set handles
set(handles.dicomcheck,'Value',options.dicomimp);
set(handles.normalize_checkbox,'Value',options.normalize.do);
if options.normalize.methodn>length(handles.normmethod,'String')
set(handles.normmethod,'Value',1);
else
set(handles.normmethod,'Value',options.normalize.methodn);
end
set(handles.normcheck,'Value',options.normalize.check);

% CT coregistration
set(handles.coregct_checkbox,'Value',options.coregct.do);
if options.coregct.methodn>length(handles.coregctmethod,'String')
set(handles.coregctmethod,'Value',1);
else
set(handles.coregctmethod,'Value',options.coregct.methodn);
end
set(handles.coregthreshs,'String',options.coregct.coregthreshs);

set(handles.coregctcheck,'Value',options.coregctcheck);



set(handles.MRCT,'Value',options.modality);

if ismember(1,options.sides)
    set(handles.right_checkbox,'Value',1);
else
        set(handles.right_checkbox,'Value',0);
end
if ismember(2,options.sides)
    set(handles.left_checkbox,'Value',1);
else
    set(handles.left_checkbox,'Value',0);
end


set(handles.doreconstruction_checkbox,'Value',options.doreconstruction);


if options.automask
    set(handles.maskwindow_txt,'String','auto')
else
    set(handles.maskwindow_txt,'String',num2str(options.maskwindow));
end
set(handles.genptatlascheck,'Value',options.atl.genpt); % generate patient specific atlases
set(handles.writeout2d_checkbox,'Value',options.d2.write);
set(handles.tdcolorscheck,'Value',options.d2.col_overlay);
set(handles.tdcontourcheck,'Value',options.d2.con_overlay);
setappdata(handles.tdcontourcolor,'color',options.d2.con_color);
set(handles.tdlabelcheck,'Value',options.d2.lab_overlay);
set(handles.bbsize,'String',num2str(options.d2.bbsize));
set(handles.manualheight_checkbox,'Value',options.manualheightcorrection);
set(handles.render_checkbox,'Value',options.d3.write);
set(handles.targetpopup,'Value',options.entrypointn);
set(handles.electrode_model_popup,'Value',options.elmodeln);
set(handles.atlassetpopup,'Value',options.atlassetn);
set(handles.exportservercheck,'Value',options.d3.autoserver);




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
storeui(handles);


% --- Executes on button press in coregct_checkbox.
function coregct_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to coregct_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of coregct_checkbox
storeui(handles);


% --- Executes on selection change in coregctmethod.
function coregctmethod_Callback(hObject, eventdata, handles)
% hObject    handle to coregctmethod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns coregctmethod contents as cell array
%        contents{get(hObject,'Value')} returns selected item from coregctmethod
methods=getappdata(gcf,'coregctmethod');
wm=get(handles.coregctmethod,'Value');

[~,~,alphasug]=eval([methods{wm},'(''probe'')']);
set(handles.coregthreshs,'String',alphasug);
storeui(handles);


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
storeui(handles);



function coregthreshs_Callback(hObject, eventdata, handles)
% hObject    handle to coregthreshs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of coregthreshs as text
%        str2double(get(hObject,'String')) returns contents of coregthreshs as a double
storeui(handles);


% --- Executes during object creation, after setting all properties.
function coregthreshs_CreateFcn(hObject, eventdata, handles)
% hObject    handle to coregthreshs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
storeui(handles);


% --- Executes on button press in coregctcheck.
function coregctcheck_Callback(hObject, eventdata, handles)
% hObject    handle to coregctcheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of coregctcheck
storeui(handles);


% --- Executes on button press in ft_checkbox.
function ft_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to ft_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ft_checkbox
storeui(handles);

% --- Executes on selection change in ftmethod.
function ftmethod_Callback(hObject, eventdata, handles)
% hObject    handle to ftmethod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns ftmethod contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ftmethod
storeui(handles);

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
storeui(handles);


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
lead_connectome;


% --- Executes on button press in specify2dwrite.
function specify2dwrite_Callback(hObject, eventdata, handles)
% hObject    handle to specify2dwrite (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ea_spec2dwrite;


% --- Executes on button press in openresultdir.
function openresultdir_Callback(hObject, eventdata, handles)
% hObject    handle to openresultdir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
    load([fileparts(which('lead')),filesep,'ea_prefs']);
    outfolder = lp.dicom.outfolder;
catch
    msgbox('Please set the working directory first!', 'Error','error');
    return;
end

if ismac
    system(['open ', outfolder]);
elseif isunix
    system(['xdg-open ', outfolder]);
elseif ispc
    system(['explorer ', outfolder]);
end

cd(outfolder);

% --- Executes on button press in viewmanual.
function viewmanual_Callback(hObject, eventdata, handles)
% hObject    handle to viewmanual (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
web('http://www.lead-dbs.org/?page_id=71', '-browser')


% --- Executes on button press in dlgroupc.
function dlgroupc_Callback(hObject, eventdata, handles)
% hObject    handle to dlgroupc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
choice = questdlg('Lead will download and install the Horn 2013 Group connectome files. This may take a while...', ...
    'Start GC Download', ...
    'OK','Cancel','OK');
earoot=[fileparts(which('lead')),filesep];
if strcmp(choice,'OK')
    disp('Downloading Horn 2013 Group Connectome');
    websave([earoot,'fibers',filesep,'gc.zip'],'http://www.lead-dbs.org/release/download.php?id=group')
    disp('Done. Installing.');
    unzip([earoot,'fibers',filesep,'gc.zip'],[earoot,'fibers',filesep]);
    disp('Done. Cleaning up.');
    delete([earoot,'fibers',filesep,'gc.zip']);
    msgbox('Download and install of the group connectome is complete.');
end


% --- Executes on selection change in vizspacepopup.
function vizspacepopup_Callback(hObject, eventdata, handles)
% hObject    handle to vizspacepopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns vizspacepopup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from vizspacepopup

if get(hObject,'Value')==2
   set(handles.writeout2d_checkbox,'Enable','off');
   set(handles.writeout2d_checkbox,'Value',0);
else
   set(handles.writeout2d_checkbox,'Enable','on');    
   %set(handles.writeout2d_checkbox,'Value',1);
end


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

outfolder=get(handles.patdir_choosebox,'String');

if strcmp(outfolder,'No Patient Selected')
    msgbox('Please set the working directory first!', 'Error','error');
    return;
end

if ismac
    system(['open ', outfolder]);
elseif isunix
    system(['xdg-open ', outfolder]);
elseif ispc
    system(['explorer ', outfolder]);
end

cd(outfolder);


% --- Executes on selection change in coregmrpopup.
function coregmrpopup_Callback(hObject, eventdata, handles)
% hObject    handle to coregmrpopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns coregmrpopup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from coregmrpopup


% --- Executes during object creation, after setting all properties.
function coregmrpopup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to coregmrpopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
