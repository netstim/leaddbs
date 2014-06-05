function varargout = eAuto(varargin)
% EAUTO M-file for eAuto.fig
%      EAUTO, by itself, creates a new EAUTO or raises the existing
%      singleton*.
%
%      H = EAUTO returns the handle to a new EAUTO or the handle to
%      the existing singleton*.
%
%      EAUTO('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in EAUTO.M with the given input arguments.
%
%      EAUTO('Property','Value',...) creates a new EAUTO or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before eAuto_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to eAuto_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help eAuto

% Last Modified by GUIDE v2.5 16-Feb-2014 18:54:20

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @eAuto_OpeningFcn, ...
                   'gui_OutputFcn',  @eAuto_OutputFcn, ...
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


% --- Executes just before eAuto is made visible.
function eAuto_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to eAuto (see VARARGIN)

% Choose default command line output for eAuto
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
ea_dispbn;


set(handles.verbositylevel_popup,'Value',3);
set(handles.normalize_checkbox,'Value',0);
%keyboard
%imshow(handles.bgimage,'bg_gui.png');

set(hObject,'Color',[1 1 1]);

set(handles.versiontxt,'String',ea_getvsn);

imshow('bg_gui.png')


% UIWAIT makes eAuto wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = eAuto_OutputFcn(hObject, eventdata, handles) 
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


%% set options

uipatdir=get(handles.patdir_choosebox,'String');

options.earoot=[fileparts(which('eAuto')),filesep];
options.root=[fileparts(uipatdir),filesep]; %'/Volumes/EspionageMounts/andreashorn/1065433271/bg/out/';
options.normalize.do=(get(handles.normalize_checkbox,'Value') == get(handles.normalize_checkbox,'Max'));
options.normalize.method=get(handles.normmethod,'Value');

options.endtolerance=str2num(get(handles.endtol_txt,'String'));
10; % how many slices to use with zero signal until end of electrode estimate.
options.sprungwert=str2num(get(handles.distance_txt,'String')); % how far electrode centroid may be (in xy axis) from last to current slice.
options.refinesteps=0; % how often to re-iterate to reconstruct trajectory. More than 2 should usually not be beneficial. Use 0 to use the direct measurement.
options.tra_stdfactor=str2num(get(handles.axis_std_txt,'String'));%0.9; % Default: 2.0 - the lower this factor, the lower the threshold (more included pixels in tra process).
options.cor_stdfactor=str2num(get(handles.cor_std_txt,'String')); % Default: 1.0 - the higher this factor, the lower the threshold (more included pixels in cor process).


% set modality (MR/CT) in options
options.modality = get(handles.MRCT,'Value');




options.verbose=get(handles.verbositylevel_popup,'Value'); % 4: Show figures but close them 3: Show all but close all figs except resultfig 2: Show all and leave figs open, 1: Show displays only, 0: Show no feedback. 
sidelog=[get(handles.left_checkbox,'Value') == get(handles.left_checkbox,'Max'),get(handles.right_checkbox,'Value') == get(handles.right_checkbox,'Max')];
sidepos=[1,2];

options.sides=sidepos(logical(sidelog)); %side=1 -> left electrode, side=2 -> right electrode. both: [1:2]

options.doreconstruction=(get(handles.doreconstruction_checkbox,'Value') == get(handles.doreconstruction_checkbox,'Max'));
if strcmp(get(handles.maskwindow_txt,'String'),'auto')
options.maskwindow=10; % initialize at 10
options.automask=1; % set automask flag
else
options.maskwindow=str2num(get(handles.maskwindow_txt,'String')); % size of the window that follows the trajectory
options.automask=0; % unset automask flag
end
options.slow=(get(handles.slowdemo_checkbox,'Value') == get(handles.slowdemo_checkbox,'Max')); % if true, there will be some pauses at critical points so that the process can be better visualized. Mainly for demonstration or debugging problems.
options.axiscontrast=(get(handles.axispopup,'Value')); % if 8: use tra only but smooth it before. % if 9: use mean of cor and tra but smooth it. % if 10: use raw tra only.
options.zheights=(get(handles.zpopup,'Value')); % if 1: use cor only, 2: use smoothed version of cor only, if 3: use mean of cor and tra, if 4: use multiplied cor * tra 5: use ^10 version of cor, 6: use ^10 version of cor.*tra 7: smoothed and then like 5. -> to determine heights of electrode contacts.
options.zresolution=10; % voxels are being parcellated into this amount of portions.
options.prolong_electrode=2*(get(handles.prolong_electrode_checkbox,'Value') == get(handles.prolong_electrode_checkbox,'Max'));
options.showatlases=(get(handles.showatlases_checkbox,'Value') == get(handles.showatlases_checkbox,'Max'));
options.showfibers=(get(handles.showfibers_checkbox,'Value') == get(handles.showfibers_checkbox,'Max'));
options.showconnectivities=(get(handles.showconnectivities_checkbox,'Value') == get(handles.showconnectivities_checkbox,'Max'));
options.writeoutimages=(get(handles.writeout2d_checkbox,'Value') == get(handles.writeout2d_checkbox,'Max'));
options.writeoutatlases=(get(handles.showatlases2d_checkbox,'Value') == get(handles.showatlases2d_checkbox,'Max'));
options.twodatlasopacity=0.5;
options.overlays={'edge',   'gp_h',     'gpi',  'morel'};
options.ocolors={'white',   'grey',    'red',   'multi'};       % use off or multi, red, blue, green, yellow, magenta, cyan, white, grey and black only.
options.oopacities=[1,    .2,         1,  	1];                 % use range 0-1
options.manualheightcorrection=(get(handles.manualheight_checkbox,'Value') == get(handles.manualheight_checkbox,'Max'));
options.render3D=(get(handles.render_checkbox,'Value') == get(handles.render_checkbox,'Max'));
options.dostimulation=(get(handles.stimcheckbox,'Value') == get(handles.stimcheckbox,'Max'));
options.numcontacts=4;
options.cg25entry=(get(handles.targetpopup,'Value') == 3);


elval = get(handles.electrode_model_popup,'Value');
string_list = get(handles.electrode_model_popup,'String');
options.elmodel=string_list{elval}; % distance of electrode contacts. 3 mm for large (VIM) electrodes, 2 mm for small (STN/GPi) electrodes.


options=ea_defaultoptions(options);

[root,thispatdir]=fileparts(uipatdir);
options.patientname=thispatdir;
ea_autocoord(options);



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


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



function z_contrast_txt_Callback(hObject, eventdata, handles)
% hObject    handle to z_contrast_txt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of z_contrast_txt as text
%        str2double(get(hObject,'String')) returns contents of z_contrast_txt as a double


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


% --- Executes on button press in render_checkbox.
function render_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to render_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of render_checkbox
if get(hObject,'Value') 
    
   set(handles.showatlases_checkbox,'Enable','on'); 
    
   set(handles.prolong_electrode_checkbox,'Enable','on'); 
   set(handles.showfibers_checkbox,'Enable','on'); 
   set(handles.showconnectivities_checkbox,'Enable','on'); 
else
    
    
   set(handles.showatlases_checkbox,'Enable','off'); 
    
   set(handles.prolong_electrode_checkbox,'Enable','off'); 
   set(handles.showfibers_checkbox,'Enable','off'); 
   set(handles.showconnectivities_checkbox,'Enable','off'); 
end

% --- Executes on button press in showatlases_checkbox.
function showatlases_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to showatlases_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of showatlases_checkbox


% --- Executes on button press in showfibers_checkbox.
function showfibers_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to showfibers_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of showfibers_checkbox


% --- Executes on button press in showconnectivities_checkbox.
function showconnectivities_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to showconnectivities_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of showconnectivities_checkbox


% --- Executes on button press in writeout2d_checkbox.
function writeout2d_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to writeout2d_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of writeout2d_checkbox
if get(hObject,'Value')
    
   set(handles.showatlases2d_checkbox,'Enable','on'); 
    
else
       set(handles.showatlases2d_checkbox,'Enable','off');
    
end

% --- Executes on button press in showatlases2d_checkbox.
function showatlases2d_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to showatlases2d_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of showatlases2d_checkbox


% --- Executes on button press in manualheight_checkbox.
function manualheight_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to manualheight_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of manualheight_checkbox


% --- Executes on button press in patdir_choosebox.
function patdir_choosebox_Callback(hObject, eventdata, handles)
% hObject    handle to patdir_choosebox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
uipatdir=uigetdir;

set(handles.patdir_choosebox,'String',uipatdir);
% check for MR-files

[~,patientname]=fileparts(uipatdir);
prefs=ea_prefs(patientname);
if exist([uipatdir,filesep,prefs.tranii],'file') && exist([uipatdir,filesep,prefs.cornii],'file')
    set(handles.statusone,'String','Normalized MR-volumes found.');
elseif ~exist([uipatdir,filesep,prefs.tranii],'file') || ~exist([uipatdir,filesep,prefs.cornii],'file')
    set(handles.statusone,'String','One or more MR-volumes missing.');
    if exist([uipatdir,filesep,prefs.tranii_unnormalized],'file') && exist([uipatdir,filesep,prefs.cornii_unnormalized],'file')
        set(handles.statusone,'String','Unnormalized MR-volumes found. Set normalize option.');
    end
end

% check for reconstructions
if exist([uipatdir,filesep,'ea_coords.fcsv'],'file') && exist([uipatdir,filesep,'ea_reconstruction.mat'],'file')
    set(handles.statustwo,'String','Fiducials and Trajectory information present in folder. Will be overwritten if "Reconstruct" is set.');
elseif exist([uipatdir,filesep,'ea_coords.fcsv'],'file') && ~exist([uipatdir,filesep,'ea_reconstruction.mat'],'file')
    set(handles.statustwo,'String','Fiducials information present in folder. Will be overwritten if "Reconstruct" is set.');
elseif ~exist([uipatdir,filesep,'ea_coords.fcsv'],'file') && exist([uipatdir,filesep,'ea_reconstruction.mat'],'file')
    set(handles.statustwo,'String','Trajectory information present in folder. Will be overwritten if "Reconstruct" is set.');
elseif ~exist([uipatdir,filesep,'ea_coords.fcsv'],'file') && ~exist([uipatdir,filesep,'ea_reconstruction.mat'],'file')
    set(handles.statustwo,'String','No reconstruction available in folder. Set "Reconstruct" to start.');
end

% --- Executes on button press in left_checkbox.
function left_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to left_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of left_checkbox


% --- Executes on button press in right_checkbox.
function right_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to right_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of right_checkbox


% --- Executes on button press in normalize_checkbox.
function normalize_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to normalize_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of normalize_checkbox



function refinesteps_txt_Callback(hObject, eventdata, handles)
% hObject    handle to refinesteps_txt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of refinesteps_txt as text
%        str2double(get(hObject,'String')) returns contents of refinesteps_txt as a double


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


% --- Executes on selection change in verbositylevel_popup.
function verbositylevel_popup_Callback(hObject, eventdata, handles)
% hObject    handle to verbositylevel_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns verbositylevel_popup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from verbositylevel_popup


% --- Executes during object creation, after setting all properties.
function verbositylevel_popup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to verbositylevel_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
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
    
    
   set(handles.endtol_txt,'Enable','on'); 
   set(handles.distance_txt,'Enable','on'); 
   set(handles.axis_std_txt,'Enable','on');    
   set(handles.axispopup,'Enable','on'); 
   set(handles.cor_std_txt,'Enable','on'); 
   set(handles.zpopup,'Enable','on'); 
   set(handles.maskwindow_txt,'Enable','on'); 
   %set(handles.manualheight_checkbox,'Enable','on'); 
   set(handles.slowdemo_checkbox,'Enable','on'); 

   
   
   
else
    
   set(handles.endtol_txt,'Enable','off'); 
   set(handles.distance_txt,'Enable','off'); 
   set(handles.axis_std_txt,'Enable','off');    
   set(handles.axispopup,'Enable','off'); 
   set(handles.cor_std_txt,'Enable','off'); 
   set(handles.zpopup,'Enable','off'); 
   set(handles.maskwindow_txt,'Enable','off'); 
   %set(handles.manualheight_checkbox,'Enable','off'); 
   set(handles.slowdemo_checkbox,'Enable','off'); 
end


% --- Executes on selection change in MRCT.
function MRCT_Callback(hObject, eventdata, handles)
% hObject    handle to MRCT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns MRCT contents as cell array
%        contents{get(hObject,'Value')} returns selected item from MRCT


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


% --- Executes on button press in cg25check.
function cg25check_Callback(hObject, eventdata, handles)
% hObject    handle to cg25check (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cg25check


% --- Executes on selection change in normmethod.
function normmethod_Callback(hObject, eventdata, handles)
% hObject    handle to normmethod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns normmethod contents as cell array
%        contents{get(hObject,'Value')} returns selected item from normmethod


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
