function varargout = lead_connectomemapper(varargin)
% LEAD_CONNECTOMEMAPPER MATLAB code for lead_connectomemapper.fig
%      LEAD_CONNECTOMEMAPPER, by itself, creates a new LEAD_CONNECTOMEMAPPER or raises the existing
%      singleton*.
%
%      H = LEAD_CONNECTOMEMAPPER returns the handle to a new LEAD_CONNECTOMEMAPPER or the handle to
%      the existing singleton*.
%
%      LEAD_CONNECTOMEMAPPER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in LEAD_CONNECTOMEMAPPER.M with the given input arguments.
%
%      LEAD_CONNECTOMEMAPPER('Property','Value',...) creates a new LEAD_CONNECTOMEMAPPER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before lead_connectomemapper_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to lead_connectomemapper_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help lead_connectomemapper

% Last Modified by GUIDE v2.5 26-Aug-2016 13:21:31

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @lead_connectomemapper_OpeningFcn, ...
                   'gui_OutputFcn',  @lead_connectomemapper_OutputFcn, ...
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


% --- Executes just before lead_connectomemapper is made visible.
function lead_connectomemapper_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to lead_connectomemapper (see VARARGIN)


earoot=ea_getearoot;
im=imread([earoot,'icons',filesep,'logo_lead_connectome.png']);
image(im);
axes(handles.logoaxes);
axis off;
axis equal;
set(handles.leadfigure,'name','Lead-Connectome Mapper','color','w');
homedir=ea_gethome;
setappdata(handles.leadfigure,'uipatdir',{homedir(1:end-1)});



ea_menu_initmenu(handles,{'cluster','prefs'});


[mdl,sf]=ea_genmodlist;
set(handles.fiberspopup,'String',mdl(sf==1));
set(handles.fmripopup,'String',mdl(sf==2));
if isempty(get(handles.fiberspopup,'String'))
    set(handles.fiberspopup,'String','No structural connectome found.');
end
if isempty(get(handles.fmripopup,'String'))
    set(handles.fmripopup,'String','No functional connectome found.');
end


set(handles.versiontxt,'String',['v',ea_getvsn('local')]);

% Choose default command line output for lead_connectomemapper
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes lead_connectomemapper wait for user response (see UIRESUME)
% uiwait(handles.leadfigure);


% --- Outputs from this function are returned to the command line.
function varargout = lead_connectomemapper_OutputFcn(hObject, eventdata, handles) 
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

p='/'; % default use root
try
p=pwd; % if possible use pwd instead (could not work if deployed)
end
try % finally use last patient parent dir if set.
earoot=ea_getearoot;
load([earoot,'ea_recentpatients.mat']);
p=fileparts(fullrpts{1});
end

[seeds,path]=uigetfile({'*.txt';'*.nii.gz';'*.nii'},'Please choose seed definition(s)...','MultiSelect','on');

if iscell(seeds)
    set(hObject,'String',['Multiple (',num2str(length(seeds)),')']);
else
    set(hObject,'String',[seeds]);
    seeds={seeds};
end

for s=1:length(seeds)
   seeds{s}=fullfile(path,seeds{s}); 
end

setappdata(hObject,'seeds',seeds);





% --- Executes on button press in run_button.
function run_button_Callback(hObject, eventdata, handles)
% hObject    handle to run_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

leadfig=handles.leadfigure;
ea_busyaction('on',leadfig,'mapper');


options=ea_handles2options(handles);
options.uipatdirs={''};
ea_run('run',options);

ea_busyaction('off',leadfig,'mapper');


% --- Executes on button press in exportcode.
function exportcode_Callback(hObject, eventdata, handles)
% hObject    handle to exportcode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

leadfig=handles.leadfigure;
ea_busyaction('on',leadfig,'mapper');

options=ea_handles2options(handles);

ea_run('export',options);

ea_busyaction('off',leadfig,'mapper');


% --- Executes on selection change in command.
function command_Callback(hObject, eventdata, handles)
% hObject    handle to command (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns command contents as cell array
%        contents{get(hObject,'Value')} returns selected item from command


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


% --- Executes on selection change in fmripopup.
function fmripopup_Callback(hObject, eventdata, handles)
% hObject    handle to fmripopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns fmripopup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from fmripopup


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
setappdata(hObject,'odir',[odir,filesep]);
set(hObject,'String',odir);


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

omask=uigetfile({'*.nii';'*.nii.gz'},'Choose output location');
setappdata(hObject,'omask',[omask]);
set(hObject,'String',omask);
