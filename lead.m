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

% Last Modified by GUIDE v2.5 02-Sep-2014 14:36:55

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
excludes={'.','..'};
asc(ismember(asc,excludes))=[];

asc{end+1}='Use none';

set(handles.atlassetpopup,'String',asc);


set(handles.normalize_checkbox,'Value',0);

set(hObject,'Color',[1 1 1]);
vnum=ea_getvsn('local');
set(handles.versiontxt,'String',['v ',num2str(vnum(1))]);






im = imread('bg_gui.png');
image(im);
axis off;
axis fill
warning('off');
set(handles.dicompanel,'BackgroundColor','none');
set(handles.gopanel,'BackgroundColor','none');
set(handles.psapanel,'BackgroundColor','none');
set(handles.normpanel,'BackgroundColor','none');
set(handles.reconpanel,'BackgroundColor','none');
set(handles.reviewpanel,'BackgroundColor','none');
set(handles.vizpanel,'BackgroundColor','none');
warning('on');

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
        cnt=cnt+1;
        end
    end
end
setappdata(gcf,'normmethod',normmethod);
set(handles.normmethod,'String',ndc);


ea_firstrun(handles);
getui(handles);




% UIWAIT makes lead wait for user response (see UIRESUME)
% uiwait(handles.figure1);


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

options=handles2options(handles);

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

if length(uipatdir)>1
    set(handles.patdir_choosebox,'String',['Multiple (',num2str(length(uipatdir)),')']);
else
set(handles.patdir_choosebox,'String',uipatdir{1});
end

% store patient directories in figure
setappdata(gcf,'uipatdir',uipatdir);
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



function refinesteps_txt_Callback(hObject, eventdata, handles)
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
    
    
   set(handles.axispopup,'Enable','on'); 
   set(handles.maskwindow_txt,'Enable','on'); 

   
   
   
else
    
   set(handles.axispopup,'Enable','off'); 
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


% --- Executes on button press in normfiberscheckbox.
function normfiberscheckbox_Callback(hObject, eventdata, handles)
% hObject    handle to normfiberscheckbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of normfiberscheckbox
storeui(handles);


% --- Executes on button press in tdcolorscheck.
function tdcolorscheck_Callback(hObject, eventdata, handles)
% hObject    handle to tdcolorscheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tdcolorscheck
storeui(handles);


% --- Executes on button press in tdcontourcheck.
function tdcontourcheck_Callback(hObject, eventdata, handles)
% hObject    handle to tdcontourcheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tdcontourcheck
storeui(handles);



function bbsize_Callback(hObject, eventdata, handles)
% hObject    handle to bbsize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of bbsize as text
%        str2double(get(hObject,'String')) returns contents of bbsize as a double
storeui(handles);


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


% --- Executes on button press in tdlabelcheck.
function tdlabelcheck_Callback(hObject, eventdata, handles)
% hObject    handle to tdlabelcheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tdlabelcheck
storeui(handles);


% --- Executes on button press in tdcontourcolor.
function tdcontourcolor_Callback(hObject, eventdata, handles)
% hObject    handle to tdcontourcolor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

tdcol=uisetcolor;
setappdata(hObject,'color',tdcol);
storeui(handles);


% --- Executes on button press in cmappushbutton.
function cmappushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to cmappushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
colormapeditor;
storeui(handles);


% --- Executes on button press in autoimprovecheck.
function autoimprovecheck_Callback(hObject, eventdata, handles)
% hObject    handle to autoimprovecheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of autoimprovecheck
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

chooseboxname=get(handles.patdir_choosebox,'String');

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

% check for MR-files
[~,patientname]=fileparts(uipatdir{1});
prefs=ea_prefs(patientname);

set(handles.statusone,'String','One or more MR-/CT-volumes missing.');
modality=zeros(2,1);
% check for unnormalized images first:
% CT
if exist([uipatdir{1},filesep,prefs.rawctnii_unnormalized],'file')
    set(handles.statusone,'String','Raw CT file found. Please coregister CT to preoperative MR first.');
    modality(2)=1;
end
if exist([uipatdir{1},filesep,prefs.ctnii_coregistered],'file')
    set(handles.statusone,'String','Coregistered CT file found. Set normalize option.');
    modality(2)=1;
end
% MR
if exist([uipatdir{1},filesep,prefs.tranii_unnormalized],'file') && exist([uipatdir{1},filesep,prefs.cornii_unnormalized],'file')
    set(handles.statusone,'String','Unnormalized MR-volumes found. Set normalize option.');
    modality(1)=1;
end
% MR & CT
if exist([uipatdir{1},filesep,prefs.tranii_unnormalized],'file') && exist([uipatdir{1},filesep,prefs.cornii_unnormalized],'file') && exist([uipatdir{1},filesep,prefs.rawctnii_unnormalized],'file')
    set(handles.statusone,'String','Unnormalized MR-volumes and raw CT volume found. Set normalize option or coregister CT.');
    modality(1:2)=1;
end
if exist([uipatdir{1},filesep,prefs.tranii_unnormalized],'file') && exist([uipatdir{1},filesep,prefs.cornii_unnormalized],'file') && exist([uipatdir{1},filesep,prefs.ctnii_coregistered],'file')
    set(handles.statusone,'String','Unnormalized MR-volumes and coregistered CT file found. Set normalize option.');
    modality(1:2)=1;
end


% Now check for normalized images:
if exist([uipatdir{1},filesep,prefs.ctnii],'file')
    try
    load([uipatdir{1},filesep,'ea_normmethod_applied.mat'])
    nmused=eval([norm_method_applied,'(''check'')']);
    set(handles.statusone,'String',['CT-volumes have been normalized with: ',nmused]);
    catch
    set(handles.statusone,'String','Normalized CT-volumes found.');
    end
    modality(2)=1;
end
if exist([uipatdir{1},filesep,prefs.tranii],'file') && exist([uipatdir{1},filesep,prefs.cornii],'file')
    try
    load([uipatdir{1},filesep,'ea_normmethod_applied.mat'])
    nmused=eval([norm_method_applied,'(''check'')']);
    set(handles.statusone,'String',['MR-volumes have been normalized with: ',nmused]);
    catch
    set(handles.statusone,'String','Normalized MR-volumes found.');
    end
    modality(1)=1;
end
if exist([uipatdir{1},filesep,prefs.tranii],'file') && exist([uipatdir{1},filesep,prefs.cornii],'file') && exist([uipatdir{1},filesep,prefs.ctnii],'file')
    set(handles.statusone,'String','Normalized MR- and CT-volumes found.');
    modality(1:2)=1;
end


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



% set modality (MR/CT) in options
options.modality = get(handles.MRCT,'Value');




options.verbose=3; % 4: Show figures but close them 3: Show all but close all figs except resultfig 2: Show all and leave figs open, 1: Show displays only, 0: Show no feedback. 
sidelog=[get(handles.right_checkbox,'Value') == get(handles.right_checkbox,'Max'),get(handles.left_checkbox,'Value') == get(handles.left_checkbox,'Max')];
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
options.autoimprove=(get(handles.autoimprovecheck,'Value') == get(handles.autoimprovecheck,'Max')); % if true, there will be some pauses at critical points so that the process can be better visualized. Mainly for demonstration or debugging problems.

options.axiscontrast=(get(handles.axispopup,'Value')); % if 8: use tra only but smooth it before. % if 9: use mean of cor and tra but smooth it. % if 10: use raw tra only.
options.zresolution=10; % voxels are being parcellated into this amount of portions.

options.atl.genpt=get(handles.genptatlascheck,'Value'); % generate patient specific atlases
options.atl.normalize=get(handles.normptatlascheck,'Value'); % normalize patient specific atlasset.
options.atl.can=get(handles.canatlcheck,'Value'); % display canonical atlases
options.atl.pt=get(handles.patatlcheck,'Value'); % display patient specific atlases


options.d2.write=(get(handles.writeout2d_checkbox,'Value') == get(handles.writeout2d_checkbox,'Max'));
options.d2.atlasopacity=0.5;
options.d2.col_overlay=get(handles.tdcolorscheck,'Value');
options.d2.con_overlay=get(handles.tdcontourcheck,'Value');
options.d2.con_color=getappdata(handles.tdcontourcolor,'color');
if isempty(options.d2.con_color)
    options.d2.con_color=[0,0,0]; % black
end

options.d2.lab_overlay=get(handles.tdlabelcheck,'Value');


options.d2.bbsize=str2double(get(handles.bbsize,'String'));



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

options.fiberthresh=1;
options.writeoutstats=1;

options.normalize_fibers=get(handles.normfiberscheckbox,'Value');
options.colormap=colormap;





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
set(handles.autoimprovecheck,'Value',options.autoimprove);
set(handles.axispopup,'Value',options.axiscontrast); % if 8: use tra only but smooth it before. % if 9: use mean of cor and tra but smooth it. % if 10: use raw tra only.
set(handles.genptatlascheck,'Value',options.atl.genpt); % generate patient specific atlases
set(handles.normptatlascheck,'Value',options.atl.normalize); % normalize patient specific atlasset.
set(handles.canatlcheck,'Value',options.atl.can); % display canonical atlases
set(handles.patatlcheck,'Value',options.atl.pt); % display patient specific atlases
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
set(handles.normfiberscheckbox,'Value',options.normalize_fibers);
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
