function varargout = ea_nbs_advanced(varargin)
% EA_NBS_ADVANCED MATLAB code for ea_nbs_advanced.fig
%      EA_NBS_ADVANCED, by itself, creates a new EA_NBS_ADVANCED or raises the existing
%      singleton*.
%
%      H = EA_NBS_ADVANCED returns the handle to a new EA_NBS_ADVANCED or the handle to
%      the existing singleton*.
%
%      EA_NBS_ADVANCED('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in EA_NBS_ADVANCED.M with the given input arguments.
%
%      EA_NBS_ADVANCED('Property','Value',...) creates a new EA_NBS_ADVANCED or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ea_nbs_advanced_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ea_nbs_advanced_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ea_nbs_advanced

% Last Modified by GUIDE v2.5 20-Jun-2016 17:21:37

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ea_nbs_advanced_OpeningFcn, ...
                   'gui_OutputFcn',  @ea_nbs_advanced_OutputFcn, ...
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


% --- Executes just before ea_nbs_advanced is made visible.
function ea_nbs_advanced_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ea_nbs_advanced (see VARARGIN)

earoot=[ea_getearoot];
setappdata(handles.nbsadvanced,'earoot',earoot);
set(handles.nbsadvanced,'name','NBS advanced settings','color','w');



% update UI:

try
    lc=load([earoot,'connectomics',filesep,'lc_options.mat']);
catch

    lc=ea_initlcopts([]);
end
if ~isfield(lc,'nbs') % compatibility with older stored userdata (<v1.4.9)
    lc=ea_initlcopts([],lc); % will merely add the nbs stuff
end
save([earoot,'connectomics',filesep,'lc_options.mat'],'-struct','lc');

nbsadv2handles(lc,handles);



% Choose default command line output for ea_nbs_advanced
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ea_nbs_advanced wait for user response (see UIRESUME)
% uiwait(handles.nbsadvanced);


% --- Outputs from this function are returned to the command line.
function varargout = ea_nbs_advanced_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function numpermutations_Callback(hObject, eventdata, handles)
% hObject    handle to numpermutations (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of numpermutations as text
%        str2double(get(hObject,'String')) returns contents of numpermutations as a double


% --- Executes during object creation, after setting all properties.
function numpermutations_CreateFcn(hObject, eventdata, handles)
% hObject    handle to numpermutations (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in component.
function component_Callback(hObject, eventdata, handles)
% hObject    handle to component (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns component contents as cell array
%        contents{get(hObject,'Value')} returns selected item from component


% --- Executes during object creation, after setting all properties.
function component_CreateFcn(hObject, eventdata, handles)
% hObject    handle to component (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in nbsmethod.
function nbsmethod_Callback(hObject, eventdata, handles)
% hObject    handle to nbsmethod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns nbsmethod contents as cell array
%        contents{get(hObject,'Value')} returns selected item from nbsmethod


% --- Executes during object creation, after setting all properties.
function nbsmethod_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nbsmethod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function alpha_Callback(hObject, eventdata, handles)
% hObject    handle to alpha (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of alpha as text
%        str2double(get(hObject,'String')) returns contents of alpha as a double


% --- Executes during object creation, after setting all properties.
function alpha_CreateFcn(hObject, eventdata, handles)
% hObject    handle to alpha (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in savebutton.
function savebutton_Callback(hObject, eventdata, handles)
% hObject    handle to savebutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
earoot=getappdata(handles.nbsadvanced,'earoot');
try
    lc=load([earoot,'connectomics',filesep,'lc_options.mat']);
catch
    lc=ea_initlcopts([]);
end
lc.nbs.adv.method= get(handles.nbsmethod,'Value');
lc.nbs.adv.compsize=get(handles.component,'Value');
lc.nbs.adv.perm=str2double(get(handles.numpermutations,'String'));
lc.nbs.adv.alpha=get(handles.alpha,'String');
lc.nbs.adv.exch=getappdata(handles.nbsadvanced,'exchange');

save([earoot,'connectomics',filesep,'lc_options.mat'],'-struct','lc');
close(handles.nbsadvanced);

% --- Executes on button press in exchangebutton.
function exchangebutton_Callback(hObject, eventdata, handles)
% hObject    handle to exchangebutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[pth,fi]=uigetfile({'*.mat','*.txt'},'Choose exchange block spec file...');
if ~fi % user pressed cancel
    set(hObject,'BackgroundColor',[0.93,0.93,0.93]);
    setappdata(handles.nbsadvanced,'exchange','');
    return
end
setappdata(handles.nbsadvanced,'exchange',[pth,filesep,fi]);
set(hObject,'BackgroundColor',[0.1;0.8;0.1]);



function handles=nbsadv2handles(lc,handles)
try set(handles.nbsmethod,'Value',lc.nbs.adv.method); end
try set(handles.component,'Value',lc.nbs.adv.compsize); end
try set(handles.numpermutations,'String',num2str(lc.nbs.adv.perm)); end
try set(handles.alpha,'String',lc.nbs.adv.alpha); end
try setappdata(handles.nbsadvanced,'exchange',lc.nbs.adv.exch); end


% --- Executes during object creation, after setting all properties.
function nbsadvanced_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nbsadvanced (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
label='https://sites.google.com/site/bctnet/comparison/nbs';
url='https://sites.google.com/site/bctnet/comparison/nbs';
position=[227, 17, 330,17];
ea_hyperlink_label(label, url, position);
