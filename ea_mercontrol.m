function varargout = ea_mercontrol(varargin)
%EA_MERCONTROL MATLAB code file for ea_mercontrol.fig
%      EA_MERCONTROL, by itself, creates a new EA_MERCONTROL or raises the existing
%      singleton*.
%
%      H = EA_MERCONTROL returns the handle to a new EA_MERCONTROL or the handle to
%      the existing singleton*.
%
%      EA_MERCONTROL('Property','Value',...) creates a new EA_MERCONTROL using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to ea_mercontrol_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      EA_MERCONTROL('CALLBACK') and EA_MERCONTROL('CALLBACK',hObject,...) call the
%      local function named CALLBACK in EA_MERCONTROL.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ea_mercontrol

% Last Modified by GUIDE v2.5 04-Apr-2017 12:39:19

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ea_mercontrol_OpeningFcn, ...
                   'gui_OutputFcn',  @ea_mercontrol_OutputFcn, ...
                   'gui_LayoutFcn',  [], ...
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


% --- Executes just before ea_mercontrol is made visible.
function ea_mercontrol_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)
clc
set(hObject,'Name','MER Control','CloseRequestFcn',@closefunciton,'KeyPressFcn',@ea_keypress)
% Choose default command line output for ea_mercontrol
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ea_mercontrol wait for user response (see UIRESUME)
% uiwait(handles.mercontrolfig);

resultfig=varargin{1};
options=varargin{2};
setappdata(hObject,'resultfig',resultfig);
setappdata(hObject,'options',options);
clearkeyboardcontrolgroup(handles,hObject,resultfig,options)

clearmertrajectories(handles,resultfig,options)
[~,trajectory,markers]=ea_load_reconstruction(options);
coords_mm=ea_resolvecoords(markers,options);
merstruct.length = options.prefs.mer.length; %default is 24mm
merstruct.offset = options.prefs.mer.offset; % default distance between mer tracts is 2mm
merstruct.colormap = [0.5,0,0;0.5,0.5,0;0,0.5,0;0.5,0,0.5;0,0.5,0.5;0,0,0.5]; %Maroon,Olive,Green,Purple,Teal,Navy

% Set mermarkers
mermarkers = getappdata(resultfig,'mermarkers');
if isempty(mermarkers)
    mermarkers = struct('side',{},'tract',{},'depth',{},'markertype',{},'session',{},'dat',{},'tag',{},'handle',{},'notes',{});
    setappdata(resultfig,'mermarkers',mermarkers)
end
% Set default implanted tract
set(handles.popupimplantedtract_left,'Value',options.prefs.mer.defaulttract+1)
set(handles.popupimplantedtract_right,'Value',options.prefs.mer.defaulttract+1)
set(handles.editimplanteddepth_left,'String','0')
set(handles.editimplanteddepth_right,'String','0')

setappdata(resultfig,'merstruct',merstruct)
merstruct.defaultmer = ea_coords2mer(coords_mm,trajectory,resultfig,options);
clear coords_mm trajectory markers

set(0,'CurrentFigure',resultfig); hold on;
color=merstruct.colormap;
for side = options.sides
    % contact_spacing = getfield(getappdata(resultfig,'elspec'),'contact_spacing');
    
    merstruct.defaultmer.central.trajectory{side} = ea_getmertrajectory(merstruct.defaultmer.central.coords_mm{side},0,merstruct.length,50);
    merstruct.defaultmer.anterior.trajectory{side} = ea_getmertrajectory(merstruct.defaultmer.anterior.coords_mm{side},0,merstruct.length,50);
    merstruct.defaultmer.posterior.trajectory{side} = ea_getmertrajectory(merstruct.defaultmer.posterior.coords_mm{side},0,merstruct.length,50);
    merstruct.defaultmer.lateral.trajectory{side} = ea_getmertrajectory(merstruct.defaultmer.lateral.coords_mm{side},0,merstruct.length,50);
    merstruct.defaultmer.medial.trajectory{side} = ea_getmertrajectory(merstruct.defaultmer.medial.coords_mm{side},0,merstruct.length,50);
    
    if side==1 % right = 1 
        merhandles.central_right = plot3(merstruct.defaultmer.central.trajectory{side}(:,1),merstruct.defaultmer.central.trajectory{side}(:,2),merstruct.defaultmer.central.trajectory{side}(:,3),'color',color(1,:),'linew',5);
        merhandles.anterior_right = plot3(merstruct.defaultmer.anterior.trajectory{side}(:,1),merstruct.defaultmer.anterior.trajectory{side}(:,2),merstruct.defaultmer.anterior.trajectory{side}(:,3),'color',color(2,:),'linew',5);
        merhandles.posterior_right = plot3(merstruct.defaultmer.posterior.trajectory{side}(:,1),merstruct.defaultmer.posterior.trajectory{side}(:,2),merstruct.defaultmer.posterior.trajectory{side}(:,3),'color',color(3,:),'linew',5);
        merhandles.lateral_right = plot3(merstruct.defaultmer.lateral.trajectory{side}(:,1),merstruct.defaultmer.lateral.trajectory{side}(:,2),merstruct.defaultmer.lateral.trajectory{side}(:,3),'color',color(4,:),'linew',5);
        merhandles.medial_right = plot3(merstruct.defaultmer.medial.trajectory{side}(:,1),merstruct.defaultmer.medial.trajectory{side}(:,2),merstruct.defaultmer.medial.trajectory{side}(:,3),'color',color(5,:),'linew',5);
    elseif side==2 % left = 2
        merhandles.central_left = plot3(merstruct.defaultmer.central.trajectory{side}(:,1),merstruct.defaultmer.central.trajectory{side}(:,2),merstruct.defaultmer.central.trajectory{side}(:,3),'color',color(1,:),'linew',5);
        merhandles.anterior_left = plot3(merstruct.defaultmer.anterior.trajectory{side}(:,1),merstruct.defaultmer.anterior.trajectory{side}(:,2),merstruct.defaultmer.anterior.trajectory{side}(:,3),'color',color(2,:),'linew',5);
        merhandles.posterior_left = plot3(merstruct.defaultmer.posterior.trajectory{side}(:,1),merstruct.defaultmer.posterior.trajectory{side}(:,2),merstruct.defaultmer.posterior.trajectory{side}(:,3),'color',color(3,:),'linew',5);
        merhandles.lateral_left = plot3(merstruct.defaultmer.lateral.trajectory{side}(:,1),merstruct.defaultmer.lateral.trajectory{side}(:,2),merstruct.defaultmer.lateral.trajectory{side}(:,3),'color',color(4,:),'linew',5);
        merhandles.medial_left = plot3(merstruct.defaultmer.medial.trajectory{side}(:,1),merstruct.defaultmer.medial.trajectory{side}(:,2),merstruct.defaultmer.medial.trajectory{side}(:,3),'color',color(5,:),'linew',5);    
    end
    
end

merstruct.currentmer=merstruct.defaultmer;
setappdata(resultfig,'merstruct',merstruct)
setappdata(resultfig,'merhandles',merhandles)
getsettogglestates(handles);



% --- Outputs from this function are returned to the command line.
function varargout = ea_mercontrol_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in togglecentral_left.
function togglecentral_left_Callback(hObject, eventdata, handles)
% hObject    handle to togglecentral_left (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglecentral_left
getsettogglestates(handles);


% --- Executes on button press in toggleanterior_left.
function toggleanterior_left_Callback(hObject, eventdata, handles)
% hObject    handle to toggleanterior_left (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of toggleanterior_left
getsettogglestates(handles);


% --- Executes on button press in toggleposterior_left.
function toggleposterior_left_Callback(hObject, eventdata, handles)
% hObject    handle to toggleposterior_left (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of toggleposterior_left
getsettogglestates(handles);


% --- Executes on button press in togglelateral_left.
function togglelateral_left_Callback(hObject, eventdata, handles)
% hObject    handle to togglelateral_left (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglelateral_left
getsettogglestates(handles);

% --- Executes on button press in togglemedial_left.
function togglemedial_left_Callback(hObject, eventdata, handles)
% hObject    handle to togglemedial_left (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglemedial_left
getsettogglestates(handles);


function posmedial_left_Callback(hObject, eventdata, handles)
% hObject    handle to posmedial_left (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of posmedial_left as text
%        str2double(get(hObject,'String')) returns contents of posmedial_left as a double
side = 2; %left
resultfig=getappdata(handles.mercontrolfig,'resultfig');
merstruct=getappdata(resultfig,'merstruct');
dist = str2double(get(hObject,'String'))-str2double(get(handles.editimplanteddepth_left,'String'));
trajectory = ea_getmertrajectory(merstruct.currentmer.medial.coords_mm{side},dist,merstruct.length,50);
ea_updatemertrajectory(handles,trajectory,'medial_left')


% --- Executes during object creation, after setting all properties.
function posmedial_left_CreateFcn(hObject, eventdata, handles)
% hObject    handle to posmedial_left (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function posposterior_right_Callback(hObject, eventdata, handles)
% hObject    handle to posposterior_right (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of posposterior_right as text
%        str2double(get(hObject,'String')) returns contents of posposterior_right as a double
side = 1; %right
resultfig=getappdata(handles.mercontrolfig,'resultfig');
merstruct=getappdata(resultfig,'merstruct');
dist = str2double(get(hObject,'String'))-str2double(get(handles.editimplanteddepth_right,'String'));
trajectory = ea_getmertrajectory(merstruct.currentmer.posterior.coords_mm{side},dist,merstruct.length,50);
ea_updatemertrajectory(handles,trajectory,'posterior_right')


% --- Executes during object creation, after setting all properties.
function posposterior_right_CreateFcn(hObject, eventdata, handles)
% hObject    handle to posposterior_right (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in togglecentral_right.
function togglecentral_right_Callback(hObject, eventdata, handles)
% hObject    handle to togglecentral_right (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglecentral_right
getsettogglestates(handles);


% --- Executes on button press in toggleposterior_right.
function toggleposterior_right_Callback(hObject, eventdata, handles)
% hObject    handle to toggleposterior_right (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of toggleposterior_right
getsettogglestates(handles);


% --- Executes on button press in toggleanterior_right.
function toggleanterior_right_Callback(hObject, eventdata, handles)
% hObject    handle to toggleanterior_right (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of toggleanterior_right
getsettogglestates(handles);


% --- Executes on button press in togglelateral_right.
function togglelateral_right_Callback(hObject, eventdata, handles)
% hObject    handle to togglelateral_right (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglelateral_right
getsettogglestates(handles);


% --- Executes on button press in togglemedial_right.
function togglemedial_right_Callback(hObject, eventdata, handles)
% hObject    handle to togglemedial_right (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglemedial_right
getsettogglestates(handles);



function poscentral_left_Callback(hObject, eventdata, handles)
% hObject    handle to poscentral_left (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of poscentral_left as text
%        str2double(get(hObject,'String')) returns contents of poscentral_left as a double
side = 2; %left
resultfig=getappdata(handles.mercontrolfig,'resultfig');
merstruct=getappdata(resultfig,'merstruct');
dist = str2double(get(hObject,'String'))-str2double(get(handles.editimplanteddepth_left,'String'));
trajectory = ea_getmertrajectory(merstruct.currentmer.central.coords_mm{side},dist,merstruct.length,50);
ea_updatemertrajectory(handles,trajectory,'central_left')


% --- Executes during object creation, after setting all properties.
function poscentral_left_CreateFcn(hObject, eventdata, handles)
% hObject    handle to poscentral_left (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function posanterior_left_Callback(hObject, eventdata, handles)
% hObject    handle to posanterior_left (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of posanterior_left as text
%        str2double(get(hObject,'String')) returns contents of posanterior_left as a double
side = 2; %left
resultfig=getappdata(handles.mercontrolfig,'resultfig');
merstruct=getappdata(resultfig,'merstruct');
dist = str2double(get(hObject,'String'))-str2double(get(handles.editimplanteddepth_left,'String'));
trajectory = ea_getmertrajectory(merstruct.currentmer.anterior.coords_mm{side},dist,merstruct.length,50);
ea_updatemertrajectory(handles,trajectory,'anterior_left')


% --- Executes during object creation, after setting all properties.
function posanterior_left_CreateFcn(hObject, eventdata, handles)
% hObject    handle to posanterior_left (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function posposterior_left_Callback(hObject, eventdata, handles)
% hObject    handle to posposterior_left (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of posposterior_left as text
%        str2double(get(hObject,'String')) returns contents of posposterior_left as a double
side = 2; %left
resultfig=getappdata(handles.mercontrolfig,'resultfig');
merstruct=getappdata(resultfig,'merstruct');
dist = str2double(get(hObject,'String'))-str2double(get(handles.editimplanteddepth_left,'String'));
trajectory = ea_getmertrajectory(merstruct.currentmer.posterior.coords_mm{side},dist,merstruct.length,50);
ea_updatemertrajectory(handles,trajectory,'posterior_left')


% --- Executes during object creation, after setting all properties.
function posposterior_left_CreateFcn(hObject, eventdata, handles)
% hObject    handle to posposterior_left (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function poslateral_left_Callback(hObject, eventdata, handles)
% hObject    handle to poslateral_left (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of poslateral_left as text
%        str2double(get(hObject,'String')) returns contents of poslateral_left as a double
side = 2; %left
resultfig=getappdata(handles.mercontrolfig,'resultfig');
merstruct=getappdata(resultfig,'merstruct');
dist = str2double(get(hObject,'String'))-str2double(get(handles.editimplanteddepth_left,'String'));
trajectory = ea_getmertrajectory(merstruct.currentmer.lateral.coords_mm{side},dist,merstruct.length,50);
ea_updatemertrajectory(handles,trajectory,'lateral_left')


% --- Executes during object creation, after setting all properties.
function poslateral_left_CreateFcn(hObject, eventdata, handles)
% hObject    handle to poslateral_left (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in keycentral_left.
function keycentral_left_Callback(hObject, eventdata, handles)
% hObject    handle to keycentral_left (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of keycentral_left


% --- Executes on button press in keyposterior_left.
function keyposterior_left_Callback(hObject, eventdata, handles)
% hObject    handle to keyposterior_left (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of keyposterior_left


% --- Executes on button press in keymedial_left.
function keymedial_left_Callback(hObject, eventdata, handles)
% hObject    handle to keymedial_left (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of keymedial_left


% --- Executes on button press in keycentral_right.
function keycentral_right_Callback(hObject, eventdata, handles)
% hObject    handle to keycentral_right (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of keycentral_right


function poscentral_right_Callback(hObject, eventdata, handles)
% hObject    handle to poscentral_right (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of poscentral_right as text
%        str2double(get(hObject,'String')) returns contents of poscentral_right as a double
side = 1; %right
resultfig=getappdata(handles.mercontrolfig,'resultfig');
merstruct=getappdata(resultfig,'merstruct');
dist = str2double(get(hObject,'String'))-str2double(get(handles.editimplanteddepth_right,'String'));
trajectory = ea_getmertrajectory(merstruct.currentmer.central.coords_mm{side},dist,merstruct.length,50);
ea_updatemertrajectory(handles,trajectory,'central_right')


% --- Executes during object creation, after setting all properties.
function poscentral_right_CreateFcn(hObject, eventdata, handles)
% hObject    handle to poscentral_right (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function posanterior_right_Callback(hObject, eventdata, handles)
% hObject    handle to posanterior_right (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of posanterior_right as text
%        str2double(get(hObject,'String')) returns contents of posanterior_right as a double
side = 1; %right
resultfig=getappdata(handles.mercontrolfig,'resultfig');
merstruct=getappdata(resultfig,'merstruct');
dist = str2double(get(hObject,'String'))-str2double(get(handles.editimplanteddepth_right,'String'));
trajectory = ea_getmertrajectory(merstruct.currentmer.anterior.coords_mm{side},dist,merstruct.length,50);
ea_updatemertrajectory(handles,trajectory,'anterior_right')


% --- Executes during object creation, after setting all properties.
function posanterior_right_CreateFcn(hObject, eventdata, handles)
% hObject    handle to posanterior_right (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function posmedial_right_Callback(hObject, eventdata, handles)
% hObject    handle to posmedial_right (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of posmedial_right as text
%        str2double(get(hObject,'String')) returns contents of posmedial_right as a double
side = 1; %right
resultfig=getappdata(handles.mercontrolfig,'resultfig');
merstruct=getappdata(resultfig,'merstruct');
dist = str2double(get(hObject,'String'))-str2double(get(handles.editimplanteddepth_right,'String'));
trajectory = ea_getmertrajectory(merstruct.currentmer.medial.coords_mm{side},dist,merstruct.length,50);
ea_updatemertrajectory(handles,trajectory,'medial_right')


% --- Executes during object creation, after setting all properties.
function posmedial_right_CreateFcn(hObject, eventdata, handles)
% hObject    handle to posmedial_right (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function poslateral_right_Callback(hObject, eventdata, handles)
% hObject    handle to poslateral_right (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of poslateral_right as text
%        str2double(get(hObject,'String')) returns contents of poslateral_right as a double
side = 1; %right
resultfig=getappdata(handles.mercontrolfig,'resultfig');
merstruct=getappdata(resultfig,'merstruct');
dist = str2double(get(hObject,'String'))-str2double(get(handles.editimplanteddepth_right,'String'));
trajectory = ea_getmertrajectory(merstruct.currentmer.lateral.coords_mm{side},dist,merstruct.length,50);
ea_updatemertrajectory(handles,trajectory,'lateral_right')


% --- Executes during object creation, after setting all properties.
function poslateral_right_CreateFcn(hObject, eventdata, handles)
% hObject    handle to poslateral_right (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in keyposterior_right.
function keyposterior_right_Callback(hObject, eventdata, handles)
% hObject    handle to keyposterior_right (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of keyposterior_right


% --- Executes on button press in keymedial_right.
function keymedial_right_Callback(hObject, eventdata, handles)
% hObject    handle to keymedial_right (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of keymedial_right


% --- Executes on selection change in popupcentral_left.
function popupcentral_left_Callback(hObject, eventdata, handles)
% hObject    handle to popupcentral_left (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupcentral_left contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupcentral_left


% --- Executes during object creation, after setting all properties.
function popupcentral_left_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupcentral_left (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupposterior_left.
function popupposterior_left_Callback(hObject, eventdata, handles)
% hObject    handle to popupposterior_left (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupposterior_left contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupposterior_left


% --- Executes during object creation, after setting all properties.
function popupposterior_left_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupposterior_left (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmedial_left.
function popupmedial_left_Callback(hObject, eventdata, handles)
% hObject    handle to popupmedial_left (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmedial_left contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmedial_left


% --- Executes during object creation, after setting all properties.
function popupmedial_left_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmedial_left (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupcentral_right.
function popupcentral_right_Callback(hObject, eventdata, handles)
% hObject    handle to popupcentral_right (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupcentral_right contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupcentral_right


% --- Executes during object creation, after setting all properties.
function popupcentral_right_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupcentral_right (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupposterior_right.
function popupposterior_right_Callback(hObject, eventdata, handles)
% hObject    handle to popupposterior_right (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupposterior_right contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupposterior_right


% --- Executes during object creation, after setting all properties.
function popupposterior_right_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupposterior_right (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmedial_right.
function popupmedial_right_Callback(hObject, eventdata, handles)
% hObject    handle to popupmedial_right (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmedial_right contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmedial_right


% --- Executes during object creation, after setting all properties.
function popupmedial_right_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmedial_right (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on selection change in popupimplantedtract_left.
function popupimplantedtract_left_Callback(hObject, eventdata, handles)
% hObject    handle to popupimplantedtract_left (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupimplantedtract_left contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupimplantedtract_left
resultfig=getappdata(handles.mercontrolfig,'resultfig');
merstruct=getappdata(resultfig,'merstruct');

side=2; %left
tractstring=lower(get(handles.popupimplantedtract_left,'String'));
tracttag=tractstring{handles.popupimplantedtract_left.Value};

merstruct = ea_updatetracts(merstruct,side,tracttag,resultfig,handles);
ea_updatemertrajectory(handles,merstruct.currentmer.central.trajectory{side},'central_left')
ea_updatemertrajectory(handles,merstruct.currentmer.anterior.trajectory{side},'anterior_left')
ea_updatemertrajectory(handles,merstruct.currentmer.posterior.trajectory{side},'posterior_left')
ea_updatemertrajectory(handles,merstruct.currentmer.lateral.trajectory{side},'lateral_left')
ea_updatemertrajectory(handles,merstruct.currentmer.medial.trajectory{side},'medial_left')
setappdata(resultfig,'merstruct',merstruct);


% --- Executes during object creation, after setting all properties.
function popupimplantedtract_left_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupimplantedtract_left (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupimplantedtract_right.
function popupimplantedtract_right_Callback(hObject, eventdata, handles)
% hObject    handle to popupimplantedtract_right (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupimplantedtract_right contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupimplantedtract_right
resultfig=getappdata(handles.mercontrolfig,'resultfig');
merstruct=getappdata(resultfig,'merstruct');

side=1; %left
tractstring=lower(get(handles.popupimplantedtract_right,'String'));
tracttag=tractstring{handles.popupimplantedtract_right.Value};

merstruct = ea_updatetracts(merstruct,side,tracttag,resultfig,handles);
ea_updatemertrajectory(handles,merstruct.currentmer.central.trajectory{side},'central_right')
ea_updatemertrajectory(handles,merstruct.currentmer.anterior.trajectory{side},'anterior_right')
ea_updatemertrajectory(handles,merstruct.currentmer.posterior.trajectory{side},'posterior_right')
ea_updatemertrajectory(handles,merstruct.currentmer.lateral.trajectory{side},'lateral_right')
ea_updatemertrajectory(handles,merstruct.currentmer.medial.trajectory{side},'medial_right')
setappdata(resultfig,'merstruct',merstruct);


% --- Executes during object creation, after setting all properties.
function popupimplantedtract_right_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupimplantedtract_right (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editimplanteddepth_left_Callback(hObject, eventdata, handles)
% hObject    handle to editimplanteddepth_left (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editimplanteddepth_left as text
%        str2double(get(hObject,'String')) returns contents of editimplanteddepth_left as a double
side = 2; %left
dist = str2double(get(hObject,'String'));
resultfig=getappdata(handles.mercontrolfig,'resultfig');
options=getappdata(handles.mercontrolfig,'options');
merhandles = getappdata(resultfig,'merhandles');
refreshresultfig(handles,resultfig,options,side)

set(handles.editimplanteddepth_left,'String',num2str(dist));
getsettogglestates(handles,side,dist);


% --- Executes during object creation, after setting all properties.
function editimplanteddepth_left_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editimplanteddepth_left (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editimplanteddepth_right_Callback(hObject, eventdata, handles)
% hObject    handle to editimplanteddepth_right (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editimplanteddepth_right as text
%        str2double(get(hObject,'String')) returns contents of editimplanteddepth_right as a double
side = 1; %right
dist = str2double(get(hObject,'String'));
resultfig=getappdata(handles.mercontrolfig,'resultfig');
options=getappdata(handles.mercontrolfig,'options');
refreshresultfig(handles,resultfig,options,side)

set(handles.editimplanteddepth_right,'String',num2str(dist));
getsettogglestates(handles,side,dist);


% --- Executes during object creation, after setting all properties.
function editimplanteddepth_right_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editimplanteddepth_right (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when selected object is changed in keyboardcontrolgroup.
function keyboardcontrolgroup_SelectionChangedFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in keyboardcontrolgroup 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% mertoggles=getappdata(getappdata(gcf,'resultfig'),'mertoggles');
keymer = get(eventdata.NewValue,'tag');
% switch get(eventdata.NewValue,'tag')
%     case 'keycentral_left'
%            
%     case 'keyposterior_left'
%         
%     case 'keymedial_left'
%         
%     case 'keycentral_right'
% 
%     case 'keyposterior_right'
%         
%     case 'keymedial_right'
%         
% end
setappdata(getappdata(gcf,'resultfig'),'keymer',keymer);


% --- Executes on selection change in popupmermarkers_right.
function popupmermarkers_right_Callback(hObject, eventdata, handles)
% hObject    handle to popupmermarkers_right (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmermarkers_right contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmermarkers_right


% --- Executes during object creation, after setting all properties.
function popupmermarkers_right_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmermarkers_right (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmermarkers_left.
function popupmermarkers_left_Callback(hObject, eventdata, handles)
% hObject    handle to popupmermarkers_left (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmermarkers_left contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmermarkers_left


% --- Executes during object creation, after setting all properties.
function popupmermarkers_left_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmermarkers_left (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in clearall.
function clearall_Callback(hObject, eventdata, handles)
% hObject    handle to clearall (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of clearall
resultfig = getappdata(handles.mercontrolfig,'resultfig');
options = getappdata(handles.mercontrolfig,'options');
clearmertrajectories(handles,resultfig,options)
clearmermarkers(handles,resultfig,options)
set(hObject,'Value',0)


% --- Executes during object creation, after setting all properties.
function clearall_CreateFcn(hObject, eventdata, handles)
% hObject    handle to clearall (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
root = fileparts(which('ea_coregmr.m'));
background = imread([root,'/icons/delete.png']);
set(hObject,'CData',background,'Value',0)


% --- Executes on button press in setdefault.
function setdefault_Callback(hObject, eventdata, handles)
% hObject    handle to setdefault (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of setdefault
resultfig = getappdata(handles.mercontrolfig,'resultfig');
options = getappdata(handles.mercontrolfig,'options');
merstruct = getappdata(resultfig,'merstruct');
setdefaultmer(handles,resultfig,merstruct,options)
getsettogglestates(handles);
set(hObject,'Value',0)


% --- Executes during object creation, after setting all properties.
function setdefault_CreateFcn(hObject, eventdata, handles)
% hObject    handle to setdefault (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
root = fileparts(which('ea_coregmr.m'));
background = imread([root,'/icons/checkmark.png']);
set(hObject,'CData',background,'Value',0)


% --- Executes on button press in undomarker.
function undomarker_Callback(hObject, eventdata, handles)
% hObject    handle to undomarker (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of undomarker
resultfig = getappdata(handles.mercontrolfig,'resultfig');
mermarkers = getappdata(resultfig,'mermarkers');
markerstring = getappdata(resultfig,'markerstring');
n = length(mermarkers);
if n~=0
    tmpmer = mermarkers(n);
    mermarkers(n)=[];
    set(tmpmer.handle,'Visible','off')
    set(tmpmer.tag.handle,'Visible','off')
    setappdata(resultfig,'tmpmer',tmpmer);
    
    if strcmpi(tmpmer.side,'right')
        markerstring.right(end) = [];
    elseif strcmpi(tmpmer.side,'left')
        markerstring.left(end) = [];
    end
    
    setappdata(resultfig,'mermarkers',mermarkers);
    setappdata(resultfig,'markerstring',markerstring)
    set(handles.popupmermarkers_right,'Visible','on','String',markerstring.right,'Value',1)
    set(handles.popupmermarkers_left,'Visible','on','String',markerstring.left,'Value',1)

end
set(hObject,'Value',0)


% --- Executes during object creation, after setting all properties.
function undomarker_CreateFcn(hObject, eventdata, handles)
% hObject    handle to undomarker (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
root = fileparts(which('ea_coregmr.m'));
background = imread([root,'/icons/undo.png']);
set(hObject,'CData',background,'Value',0)


% --- Executes on button press in redomarker.
function redomarker_Callback(hObject, eventdata, handles)
% hObject    handle to redomarker (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of redomarker
resultfig = getappdata(handles.mercontrolfig,'resultfig');
options = getappdata(handles.mercontrolfig,'options');
mermarkers = getappdata(resultfig,'mermarkers');
markerstring = getappdata(resultfig,'markerstring');
tmpmer = getappdata(resultfig,'tmpmer');
n = length(mermarkers);
set(0,'CurrentFigure',resultfig);
if ~isempty(tmpmer)
    mermarkers(n+1)=tmpmer;
    % display marker
    tmp = tmpmer.handle;
    mermarkers(n+1).handle  = surf(tmp.XData,tmp.YData,tmp.ZData,'FaceLighting','gouraud',...
    'FaceColor',tmp.FaceColor,'EdgeColor',tmp.EdgeColor,'FaceAlpha',tmp.FaceAlpha,'Visible','on');
    clear tmp
    
    % display tag
    tmp = tmpmer.tag.handle;
    mermarkers(n+1).tag.handle = text(tmp.Position(1),tmp.Position(2),tmp.Position(3),...
    tmp.String,'Color',tmp.Color,'HorizontalAlignment',tmp.HorizontalAlignment,'Visible','on');
    clear tmp
    
    keymer = tmpmer.name;
    if strcmp(keymer(strfind(keymer,'_')+1:end),'right')
        markerstring.right{end+1} = sprintf('%0.0f. %s',n+1,tmpmer.tag.string);
    elseif strcmp(keymer(strfind(keymer,'_')+1:end),'left')
        markerstring.left{end+1} = sprintf('%0.0f. %s',n+1,tmpmer.tag.string);
    end
    
    tmp = []; %rmfield(getappdata(resultfig),'tmpmer');
    setappdata(resultfig,'tmpmer',tmp);
    setappdata(resultfig,'mermarkers',mermarkers);
    setappdata(resultfig,'markerstring',markerstring)
    set(handles.popupmermarkers_right,'Visible','on','String',markerstring.right,'Value',1)
    set(handles.popupmermarkers_left,'Visible','on','String',markerstring.left,'Value',1)

end
set(hObject,'Value',0)


% --- Executes during object creation, after setting all properties.
function redomarker_CreateFcn(hObject, eventdata, handles)
% hObject    handle to redomarker (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
root = fileparts(which('ea_coregmr.m'));
background = imread([root,'/icons/redo.png']);
set(hObject,'CData',background,'Value',0)



% --- Executes on button press in togglemarkertags.
function togglemarkertags_Callback(hObject, eventdata, handles)
% hObject    handle to togglemarkertags (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglemarkertags
resultfig = getappdata(handles.mercontrolfig,'resultfig');
options = getappdata(handles.mercontrolfig,'options');
mermarkers = getappdata(resultfig,'mermarkers');
if get(hObject,'Value')
    for i = 1:length(mermarkers)
        try set(mermarkers(i).tag.handle,'Visible','on'); end
    end
    % set(handles.togglemarkertags,'Value',1)
else
    for i = 1:length(mermarkers)
        try set(mermarkers(i).tag.handle,'Visible','off'); end
    end
    % set(handles.togglemarkertags,'Value',0)
end


% --- Executes during object creation, after setting all properties.
function togglemarkertags_CreateFcn(hObject, eventdata, handles)
% hObject    handle to togglemarkertags (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
root = fileparts(which('ea_coregmr.m'));
background = imread([root,'/icons/text.png']);
set(hObject,'CData',background,'Value',0)


% --- Executes on button press in importmarkers.
function importmarkers_Callback(hObject, eventdata, handles)
% hObject    handle to importmarkers (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of importmarkers
resultfig = getappdata(handles.mercontrolfig,'resultfig');
options = getappdata(handles.mercontrolfig,'options');
merhandles = getappdata(resultfig,'merhandles');
if isvalid(merhandles.central_right)
try 
    load([options.uipatdirs{1},filesep,'ea_mermarkers.mat'],'mermarkers','mertoggles');
catch
    [files] = fdir(options.uipatdirs{1},'.mat');
    if ~isempty(cell2mat(strfind({files(:).name},'ea_mermarkers')))
        %[s,b] = listdlg('PromptString','Select a file:','SelectionMode','single',...
        %   'ListString',{files(~cellfun(@isempty,strfind({files(:).name},'ea_mermarkers'))).name})
    end
    error('could not find ea_mermarkers.mat')
end
%assignin(workspace,'mermarkers',mermarkers)
getsettogglestates(handles,mertoggles);
ea_updatemermarkers(mermarkers,handles,resultfig,options)
else
    fprintf(2,'Error using importmarkers_Callback: \nNo MER trajectories available\n\n')
end
set(hObject,'Value',0)


% --- Executes during object creation, after setting all properties.
function importmarkers_CreateFcn(hObject, eventdata, handles)
% hObject    handle to importmarkers (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
root = fileparts(which('ea_coregmr.m'));
background = imread([root,'/icons/import.png']);
set(hObject,'CData',background,'Value',0)


% --- Executes on button press in exportmarkers.
function exportmarkers_Callback(hObject, eventdata, handles)
% hObject    handle to exportmarkers (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of exportmarkers
resultfig = getappdata(handles.mercontrolfig,'resultfig');
options = getappdata(handles.mercontrolfig,'options');
mermarkers = getappdata(resultfig,'mermarkers');
mertoggles = getsettogglestates(handles);
if exist([options.uipatdirs{1},filesep,'ea_mermarkers.mat'],'file')
    overwrite = ea_questdlg({['ea_mermarkers.mat found in ' options.patientname ' directory.'];'Would you like to overwrite this file?'},options.patientname);
end
if ~exist([options.uipatdirs{1},filesep,'ea_mermarkers.mat'],'file') || strcmp(overwrite,'Yes')
    disp(['Saving: ' options.uipatdirs{1},filesep,'ea_mermarkers.mat'])
    save([options.uipatdirs{1},filesep,'ea_mermarkers.mat'],'mermarkers','mertoggles')
    disp('DONE')
end    
set(hObject,'Value',0)


% --- Executes during object creation, after setting all properties.
function exportmarkers_CreateFcn(hObject, eventdata, handles)
% hObject    handle to exportmarkers (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
root = fileparts(which('ea_coregmr.m'));
background = imread([root,'/icons/export.png']);
set(hObject,'CData',background,'Value',0)




% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
function setdefaultmer(handles,resultfig,merstruct,options)
clearmertrajectories(handles,resultfig,options)
color=merstruct.colormap;
set(0,'CurrentFigure',resultfig)
for side = options.sides
    if side==1
        merhandles.central_right = plot3(merstruct.defaultmer.central.trajectory{side}(:,1),merstruct.defaultmer.central.trajectory{side}(:,2),merstruct.defaultmer.central.trajectory{side}(:,3),'color',color(1,:),'linew',5);
        merhandles.anterior_right = plot3(merstruct.defaultmer.anterior.trajectory{side}(:,1),merstruct.defaultmer.anterior.trajectory{side}(:,2),merstruct.defaultmer.anterior.trajectory{side}(:,3),'color',color(2,:),'linew',5);
        merhandles.posterior_right = plot3(merstruct.defaultmer.posterior.trajectory{side}(:,1),merstruct.defaultmer.posterior.trajectory{side}(:,2),merstruct.defaultmer.posterior.trajectory{side}(:,3),'color',color(3,:),'linew',5);
        merhandles.lateral_right = plot3(merstruct.defaultmer.lateral.trajectory{side}(:,1),merstruct.defaultmer.lateral.trajectory{side}(:,2),merstruct.defaultmer.lateral.trajectory{side}(:,3),'color',color(4,:),'linew',5);
        merhandles.medial_right = plot3(merstruct.defaultmer.medial.trajectory{side}(:,1),merstruct.defaultmer.medial.trajectory{side}(:,2),merstruct.defaultmer.medial.trajectory{side}(:,3),'color',color(5,:),'linew',5);
        setappdata(resultfig,'merhandles',merhandles)

    elseif side==2
        merhandles.central_left = plot3(merstruct.defaultmer.central.trajectory{side}(:,1),merstruct.defaultmer.central.trajectory{side}(:,2),merstruct.defaultmer.central.trajectory{side}(:,3),'color',color(1,:),'linew',5);
        merhandles.anterior_left = plot3(merstruct.defaultmer.anterior.trajectory{side}(:,1),merstruct.defaultmer.anterior.trajectory{side}(:,2),merstruct.defaultmer.anterior.trajectory{side}(:,3),'color',color(2,:),'linew',5);
        merhandles.posterior_left = plot3(merstruct.defaultmer.posterior.trajectory{side}(:,1),merstruct.defaultmer.posterior.trajectory{side}(:,2),merstruct.defaultmer.posterior.trajectory{side}(:,3),'color',color(3,:),'linew',5);
        merhandles.lateral_left = plot3(merstruct.defaultmer.lateral.trajectory{side}(:,1),merstruct.defaultmer.lateral.trajectory{side}(:,2),merstruct.defaultmer.lateral.trajectory{side}(:,3),'color',color(4,:),'linew',5);
        merhandles.medial_left = plot3(merstruct.defaultmer.medial.trajectory{side}(:,1),merstruct.defaultmer.medial.trajectory{side}(:,2),merstruct.defaultmer.medial.trajectory{side}(:,3),'color',color(5,:),'linew',5);    
        setappdata(resultfig,'merhandles',merhandles)
    end
end


function clearmermarkers(handles,resultfig,options)
% clear mermarkers and reset structure
mermarkers = getappdata(resultfig,'mermarkers');
for i = 1:length(mermarkers)
try delete(mermarkers(i).handle); end
try delete(mermarkers(i).tag.handle); end
end
mermarkers = struct('side',{},'tract',{},'depth',{},'markertype',{},'session',{},'dat',{},'tag',{},'handle',{},'notes',{});
setappdata(resultfig,'mermarkers',mermarkers)
set(handles.popupimplantedtract_left,'Value',options.prefs.mer.defaulttract+1);
set(handles.popupimplantedtract_right,'Value',options.prefs.mer.defaulttract+1);
set(handles.editimplanteddepth_left,'String','0');
set(handles.editimplanteddepth_right,'String','0');
set(handles.popupmermarkers_left,'Visible','off','String','','Value',1)
set(handles.popupmermarkers_right,'Visible','off','String','','Value',1)


function clearmertrajectories(handles,resultfig,options,side)
if nargin<4
    side=[1 2];
end

merhandles = getappdata(resultfig,'merhandles');
for i = side
    if i==1
    try; delete(merhandles.central_right); end
    try; delete(merhandles.anterior_right); end
    try; delete(merhandles.posterior_right); end
    try; delete(merhandles.lateral_right); end
    try; delete(merhandles.medial_right); end
    
    poscentral_right = getfield(getappdata(handles.mercontrolfig,'UsedByGUIData_m'),'poscentral_right');
    posanterior_right = getfield(getappdata(handles.mercontrolfig,'UsedByGUIData_m'),'posanterior_right');
    posposterior_right = getfield(getappdata(handles.mercontrolfig,'UsedByGUIData_m'),'posposterior_right');
    poslateral_right = getfield(getappdata(handles.mercontrolfig,'UsedByGUIData_m'),'poslateral_right');
    posmedial_right = getfield(getappdata(handles.mercontrolfig,'UsedByGUIData_m'),'posmedial_right');
    
    set(poscentral_right,'String','0');
    set(posanterior_right,'String','0');
    set(posposterior_right,'String','0');
    set(poslateral_right,'String','0');
    set(posmedial_right,'String','0');

    elseif i==2
        try; delete(merhandles.central_left); end
        try; delete(merhandles.anterior_left); end
        try; delete(merhandles.posterior_left); end
        try; delete(merhandles.lateral_left); end
        try; delete(merhandles.medial_left); end

        poscentral_left = getfield(getappdata(handles.mercontrolfig,'UsedByGUIData_m'),'poscentral_left');
        posanterior_left = getfield(getappdata(handles.mercontrolfig,'UsedByGUIData_m'),'posanterior_left');
        posposterior_left = getfield(getappdata(handles.mercontrolfig,'UsedByGUIData_m'),'posposterior_left');
        poslateral_left = getfield(getappdata(handles.mercontrolfig,'UsedByGUIData_m'),'poslateral_left');
        posmedial_left = getfield(getappdata(handles.mercontrolfig,'UsedByGUIData_m'),'posmedial_left');
        set(poscentral_left,'String','0');
        set(posanterior_left,'String','0');
        set(posposterior_left,'String','0');
        set(poslateral_left,'String','0');
        set(posmedial_left,'String','0');
    end
end
    setappdata(resultfig,'merhandles',merhandles)


function closefunciton(hObject,event)
try
resultfig = getappdata(hObject,'resultfig');
options = getappdata(hObject,'options');
handles = getappdata(hObject,'UsedByGUIData_m');
clearmertrajectories(handles,resultfig,options)
clearmermarkers(handles,resultfig,options)
end
delete(hObject)


function ea_updatemertrajectory(handles,trajectory,tag)
resultfig=getappdata(handles.mercontrolfig,'resultfig');
% Update position in resultfig
% XData = get(getappdata(resultfig,tag),'XData');
% YData = get(getappdata(resultfig,tag),'YData');
% ZData = get(getappdata(resultfig,tag),'ZData');
h = getfield(getappdata(resultfig,'merhandles'),tag);
set(h,'XData',trajectory(:,1)')
set(h,'YData',trajectory(:,2)')
set(h,'ZData',trajectory(:,3)')
setappdata(resultfig,tag,h)
set(handles.(['key' tag]),'Value',1)
setappdata(resultfig,'keymer',['key' tag])


function ea_updatemermarkers(mermarkers,handles,resultfig,options)
clearmermarkers(handles,resultfig,options)
set(0,'CurrentFigure',resultfig)
markerstring.right = {'none selected...'};
markerstring.left = {'none selected...'};
for n = 1:length(mermarkers)
    tmp = mermarkers(n).handle;
    mermarkers(n).handle = surf(tmp.XData,tmp.YData,tmp.ZData,...
        'FaceColor',tmp.FaceColor,'EdgeColor',tmp.EdgeColor,'FaceAlpha',tmp.FaceAlpha); 
    clear tmp
    tmp = mermarkers(n).tag.handle;
    mermarkers(n).tag.handle = text(tmp.Position(1),tmp.Position(2),tmp.Position(3),...
        tmp.String,'Color',tmp.Color,'HorizontalAlignment',tmp.HorizontalAlignment,'Visible','off');

    if strcmp(mermarkers(n).side,'right')
        markerstring.right{end+1} = sprintf('%0.0f. %s',n,mermarkers(n).tag.string);
    elseif strcmp(mermarkers(n).side,'left')
        markerstring.left{end+1} = sprintf('%0.0f. %s',n,mermarkers(n).tag.string);
    end
    
end

% set implanted tract data
try
    Lidx = find(~cellfun(@isempty,strfind({mermarkers.side},'left'))==1);
    set(handles.popupimplantedtract_left,'Value',find(~cellfun(@isempty,strfind(handles.popupimplantedtract_left.String,mermarkers(Lidx(1)).dat.implantedtract))==1))
    set(handles.editimplanteddepth_left,'String',num2str(mermarkers(Lidx(1)).dat.leaddepth))
    getsettogglestates(handles,2,mermarkers(Lidx(1)).dat.leaddepth);
end

try
    Ridx = ~cellfun(@isempty,strfind({mermarkers.side},'right'));
    set(handles.popupimplantedtract_right,'Value',find(~cellfun(@isempty,strfind(handles.popupimplantedtract_left.String,mermarkers(Ridx(1)).dat.implantedtract))==1))
    set(handles.editimplanteddepth_left,'String',num2str(mermarkers(Ridx(1)).dat.leaddepth))
    getsettogglestates(handles,1,mermarkers(Ridx(1)).dat.leaddepth);
end

setappdata(resultfig,'mermarkers',mermarkers)
setappdata(resultfig,'markerstring',markerstring)
%set(handles.popupmermarkers_right,'Visible','on','String',markerstring.right,'Value',1)
%set(handles.popupmermarkers_left,'Visible','on','String',markerstring.left,'Value',1)
set(handles.togglemarkertags,'Visible','on','Value',0)


function mertoggles = getsettogglestates(handles,varargin)

if size(varargin,2)==1
    mertoggles = varargin{1};
elseif size(varargin,2)==2
    side = varargin{1};
    dist = varargin{2};
end

resultfig=getappdata(handles.mercontrolfig,'resultfig');
merhandles = getappdata(resultfig,'merhandles');

% reset states based on gui:
if ~exist('mertoggles','var')
mertoggles.keycontrol=[get(handles.keycentral_left,'Value'),get(handles.keyanterior_left,'Value'),...
    get(handles.keyposterior_left,'Value'),get(handles.keylateral_left,'Value'),get(handles.keymedial_left,'Value');...
    get(handles.keycentral_right,'Value'),get(handles.keyanterior_right,'Value'),get(handles.keyposterior_right,'Value'),...
    get(handles.keylateral_right,'Value'),get(handles.keymedial_right,'Value')];
mertoggles.togglestates=[get(handles.togglecentral_left,'Value'),get(handles.toggleanterior_left,'Value'),...
    get(handles.toggleposterior_left,'Value'),get(handles.togglelateral_left,'Value'),get(handles.togglemedial_left,'Value');...
    get(handles.togglecentral_right,'Value'),get(handles.toggleanterior_right,'Value'),get(handles.toggleposterior_right,'Value'),...
    get(handles.togglelateral_right,'Value'),get(handles.togglemedial_right,'Value')];
else
set(handles.keycentral_left,'Value',mertoggles.keycontrol(1,1))
    set(handles.keyanterior_left,'Value',mertoggles.keycontrol(1,2))
    set(handles.keyposterior_left,'Value',mertoggles.keycontrol(1,3))
    set(handles.keylateral_left,'Value',mertoggles.keycontrol(1,4))
    set(handles.keymedial_left,'Value',mertoggles.keycontrol(1,5))
set(handles.keycentral_right,'Value',mertoggles.keycontrol(2,1))
    set(handles.keyanterior_right,'Value',mertoggles.keycontrol(2,2))
    set(handles.keyposterior_right,'Value',mertoggles.keycontrol(2,3))
    set(handles.keylateral_right,'Value',mertoggles.keycontrol(2,4))
    set(handles.keymedial_right,'Value',mertoggles.keycontrol(2,5))
    
set(handles.togglecentral_left,'Value',mertoggles.togglestates(1,1))
    set(handles.toggleanterior_left,'Value',mertoggles.togglestates(1,2))
    set(handles.toggleposterior_left,'Value',mertoggles.togglestates(1,3))
    set(handles.togglelateral_left,'Value',mertoggles.togglestates(1,4))
    set(handles.togglemedial_left,'Value',mertoggles.togglestates(1,5))
set(handles.togglecentral_right,'Value',mertoggles.togglestates(2,1))
    set(handles.toggleanterior_right,'Value',mertoggles.togglestates(2,2))
    set(handles.toggleposterior_right,'Value',mertoggles.togglestates(2,3))
    set(handles.togglelateral_right,'Value',mertoggles.togglestates(2,4))
    set(handles.togglemedial_right,'Value',mertoggles.togglestates(2,5))

end
setappdata(getappdata(handles.mercontrolfig,'resultfig'),'mertoggles',mertoggles); % also store toggle data in resultfig.

% set togglestates 
% toggles right
if get(handles.togglecentral_right,'Value')
    set(merhandles.central_right,'Visible','on')
elseif ~get(handles.togglecentral_right,'Value')
    set(merhandles.central_right,'Visible','off')
end
if get(handles.toggleanterior_right,'Value')
    set(merhandles.anterior_right,'Visible','on')
elseif ~get(handles.toggleanterior_right,'Value')
    set(merhandles.anterior_right,'Visible','off')
end
if get(handles.toggleposterior_right,'Value')
    set(merhandles.posterior_right,'Visible','on')
elseif ~get(handles.toggleposterior_right,'Value')
    set(merhandles.posterior_right,'Visible','off')
end
if get(handles.togglelateral_right,'Value')
    set(merhandles.lateral_right,'Visible','on')
elseif ~get(handles.togglelateral_right,'Value')
    set(merhandles.lateral_right,'Visible','off')
end
if get(handles.togglemedial_right,'Value')
    set(merhandles.medial_right,'Visible','on')
elseif ~get(handles.togglemedial_right,'Value')
    set(merhandles.medial_right,'Visible','off')
end
% toggles left
if get(handles.togglecentral_left,'Value')
    set(merhandles.central_left,'Visible','on')
elseif ~get(handles.togglecentral_left,'Value')
    set(merhandles.central_left,'Visible','off')
end
if get(handles.toggleanterior_left,'Value')
    set(merhandles.anterior_left,'Visible','on')
elseif ~get(handles.toggleanterior_left,'Value')
    set(merhandles.anterior_left,'Visible','off')
end
if get(handles.toggleposterior_left,'Value')
    set(merhandles.posterior_left,'Visible','on')
elseif ~get(handles.toggleposterior_left,'Value')
    set(merhandles.posterior_left,'Visible','off')
end
if get(handles.togglelateral_left,'Value')
    set(merhandles.lateral_left,'Visible','on')
elseif ~get(handles.togglelateral_left,'Value')
    set(merhandles.lateral_left,'Visible','off')
end
if get(handles.togglemedial_left,'Value')
    set(merhandles.medial_left,'Visible','on')
elseif ~get(handles.togglemedial_left,'Value')
    set(merhandles.medial_left,'Visible','off')
end

% set pos handles
if exist('dist','var')
    if side==2
    set(handles.poscentral_left,'String',num2str(dist));
    set(handles.posanterior_left,'String',num2str(dist));
    set(handles.posposterior_left,'String',num2str(dist));
    set(handles.poslateral_left,'String',num2str(dist));
    set(handles.posmedial_left,'String',num2str(dist));
    elseif side==1
    set(handles.poscentral_right,'String',num2str(dist));
    set(handles.posanterior_right,'String',num2str(dist));
    set(handles.posposterior_right,'String',num2str(dist));
    set(handles.poslateral_right,'String',num2str(dist));
    set(handles.posmedial_right,'String',num2str(dist));
    end
end
setappdata(resultfig,'merhandles',merhandles)



function outputtrajectory = ea_getmertrajectory(coords_mm,dist,length,n)
if size(coords_mm,1)<2
    error('Must input a vector')
end
dxyz = sqrt((diff(coords_mm(1:2,1))^2)+(diff(coords_mm(1:2,2))^2)+diff(coords_mm(1:2,3))^2);
slope = mean(diff(coords_mm))/dxyz;
startpoint = coords_mm(1,:)+slope.*dist;
%     normtrajvector=slope/norm(slope);
%     orth=null(normtrajvector);
%
%     startpoint.x = trajin(1,:)+orth(:,1)';
%     startpoint.y = trajin(1,:)+orth(:,2)';
%     startpoint=slope(1,:)-(2*(trajin(1,:)-slope(1,:)))

outputtrajectory(:,1) = linspace(startpoint(1,1),startpoint(1,1)+slope(1)*length,n);
outputtrajectory(:,2) = linspace(startpoint(1,2),startpoint(1,2)+slope(2)*length,n);
outputtrajectory(:,3) = linspace(startpoint(1,3),startpoint(1,3)+slope(3)*length,n);


function merstruct = ea_updatetracts(merstruct,side,tracttag,resultfig,handles)
        
offset = getfield(getappdata(resultfig,'merstruct'),'offset');
% elstruct = 
% x-axis --> negative = medial, positive = lateral
% y-axis --> negative = posterior, positive = anterior
% z-axis --> negative = inferior, positive = superior
coords_mm = merstruct.defaultmer.central.coords_mm{side};
trajectory = merstruct.defaultmer.central.trajectory{side};
switch tracttag 
    case 'central'
        merstruct.currentmer.central.coords_mm{side} = coords_mm;
        merstruct.currentmer.central.trajectory{side} = trajectory;
        merstruct.currentmer.anterior.trajectory{side} = [coords_mm(:,1),coords_mm(:,2)+offset,coords_mm(:,3)]; %2mm anterior
        merstruct.currentmer.anterior.trajectory{side} = [trajectory(:,1),trajectory(:,2)+offset,trajectory(:,3)]; %2mm anterior
        merstruct.currentmer.posterior.coords_mm{side} = [coords_mm(:,1),coords_mm(:,2)-offset,coords_mm(:,3)]; %2mm posterior
        merstruct.currentmer.posterior.trajectory{side} = [trajectory(:,1),trajectory(:,2)-offset,trajectory(:,3)]; %2mm posterior
        if side==1 %right
            merstruct.currentmer.lateral.coords_mm{side} = [coords_mm(:,1)+offset,coords_mm(:,2),coords_mm(:,3)]; %2mm lateral (right is positive)
            merstruct.currentmer.lateral.trajectory{side} = [trajectory(:,1)+offset,trajectory(:,2),trajectory(:,3)]; %2mm lateral (right is positive)
            merstruct.currentmer.medial.coords_mm{side} = [coords_mm(:,1)-offset,coords_mm(:,2),coords_mm(:,3)]; %2mm medial (right is positive)
            merstruct.currentmer.medial.trajectory{side} = [trajectory(:,1)-offset,trajectory(:,2),trajectory(:,3)]; %2mm medial (right is positive)
        elseif side==2 %left
            merstruct.currentmer.lateral.coords_mm{side} = [coords_mm(:,1)-offset,coords_mm(:,2),coords_mm(:,3)]; %2mm medial (left is negative)
            merstruct.currentmer.lateral.trajectory{side} = [trajectory(:,1)-offset,trajectory(:,2),trajectory(:,3)]; %2mm medial (left is negative)
            merstruct.currentmer.medial.coords_mm{side} = [coords_mm(:,1)+offset,coords_mm(:,2),coords_mm(:,3)]; %2mm medial (left is negative)
            merstruct.currentmer.medial.trajectory{side} = [trajectory(:,1)+offset,trajectory(:,2),trajectory(:,3)]; %2mm medial (left is negative)
        end
    case 'anterior'
        merstruct.currentmer.anterior.coords_mm{side} = coords_mm;
        merstruct.currentmer.anterior.trajectory{side} = trajectory;
        merstruct.currentmer.central.coords_mm{side} = [coords_mm(:,1),coords_mm(:,2)-offset,coords_mm(:,3)]; %2mm anterior
        merstruct.currentmer.central.trajectory{side} = [trajectory(:,1),trajectory(:,2)-offset,trajectory(:,3)]; %2mm anterior
        merstruct.currentmer.posterior.coords_mm{side} = [coords_mm(:,1),coords_mm(:,2)-offset*2,coords_mm(:,3)]; %4mm anterior
        merstruct.currentmer.posterior.trajectory{side} = [trajectory(:,1),trajectory(:,2)-offset*2,trajectory(:,3)]; %4mm anterior
         if side==1 %right
            merstruct.currentmer.lateral.coords_mm{side} =  [coords_mm(:,1)+offset,coords_mm(:,2)-offset,coords_mm(:,3)]; %2mm medial (right is positive)
            merstruct.currentmer.lateral.trajectory{side} = [trajectory(:,1)+offset,trajectory(:,2)-offset,trajectory(:,3)]; %2mm medial (right is positive)
            merstruct.currentmer.medial.coords_mm{side} = [coords_mm(:,1)-offset,coords_mm(:,2)-offset,coords_mm(:,3)]; %2mm medial (right is positive)
            merstruct.currentmer.medial.trajectory{side} = [trajectory(:,1)-offset,trajectory(:,2)-offset,trajectory(:,3)]; %2mm medial (right is positive)
        elseif side==2 %left
            merstruct.currentmer.lateral.coords_mm{side} = [coords_mm(:,1)-offset,coords_mm(:,2)-offset,coords_mm(:,3)]; %2mm medial (left is negative)
            merstruct.currentmer.lateral.trajectory{side} = [trajectory(:,1)-offset,trajectory(:,2)-offset,trajectory(:,3)]; %2mm medial (left is negative)
            merstruct.currentmer.medial.coords_mm{side} = [coords_mm(:,1)+offset,coords_mm(:,2)-offset,coords_mm(:,3)]; %2mm medial (left is negative)
            merstruct.currentmer.medial.trajectory{side} = [trajectory(:,1)+offset,trajectory(:,2)-offset,trajectory(:,3)]; %2mm medial (left is negative)
         end
    case 'posterior'
        merstruct.currentmer.posterior.coords_mm{side} = coords_mm;
        merstruct.currentmer.posterior.trajectory{side} = trajectory;
        merstruct.currentmer.central.coords_mm{side} = [coords_mm(:,1),coords_mm(:,2)+offset,coords_mm(:,3)]; %2mm anterior
        merstruct.currentmer.central.trajectory{side} = [trajectory(:,1),trajectory(:,2)+offset,trajectory(:,3)]; %2mm anterior
        merstruct.currentmer.anterior.coords_mm{side} = [coords_mm(:,1),coords_mm(:,2)+offset*2,coords_mm(:,3)]; %4mm anterior
        merstruct.currentmer.anterior.trajectory{side} = [trajectory(:,1),trajectory(:,2)+offset*2,trajectory(:,3)]; %4mm anterior
         if side==1 %right
            merstruct.currentmer.lateral.coords_mm{side} =  [coords_mm(:,1)+offset,coords_mm(:,2)+offset,coords_mm(:,3)]; %2mm medial (right is positive)
            merstruct.currentmer.lateral.trajectory{side} = [trajectory(:,1)+offset,trajectory(:,2)+offset,trajectory(:,3)]; %2mm medial (right is positive)
            merstruct.currentmer.medial.coords_mm{side} = [coords_mm(:,1)-offset,coords_mm(:,2)+offset,coords_mm(:,3)]; %2mm medial (right is positive)
            merstruct.currentmer.medial.trajectory{side} = [trajectory(:,1)-offset,trajectory(:,2)+offset,trajectory(:,3)]; %2mm medial (right is positive)
        elseif side==2 %left
            merstruct.currentmer.lateral.coords_mm{side} = [coords_mm(:,1)-offset,coords_mm(:,2)+offset,coords_mm(:,3)]; %2mm medial (left is negative)
            merstruct.currentmer.lateral.trajectory{side} = [trajectory(:,1)-offset,trajectory(:,2)+offset,trajectory(:,3)]; %2mm medial (left is negative)
            merstruct.currentmer.medial.coords_mm{side} = [coords_mm(:,1)+offset,coords_mm(:,2)+offset,coords_mm(:,3)]; %2mm medial (left is negative)
            merstruct.currentmer.medial.trajectory{side} = [trajectory(:,1)+offset,trajectory(:,2)+offset,trajectory(:,3)]; %2mm medial (left is negative)
         end
    case 'lateral'
        merstruct.currentmer.lateral.coords_mm{side} = coords_mm;
        merstruct.currentmer.lateral.trajectory{side} = trajectory;
        if side==1 %right
            merstruct.currentmer.central.coords_mm{side} = [coords_mm(:,1)-offset,coords_mm(:,2),coords_mm(:,3)];
            merstruct.currentmer.central.trajectory{side} = [trajectory(:,1)-offset,trajectory(:,2),trajectory(:,3)];
            merstruct.currentmer.anterior.coords_mm{side} = [coords_mm(:,1)-offset,coords_mm(:,2)+offset,coords_mm(:,3)];
            merstruct.currentmer.anterior.trajectory{side} = [trajectory(:,1)-offset,trajectory(:,2)+offset,trajectory(:,3)];
            merstruct.currentmer.posterior.coords_mm{side} = [coords_mm(:,1)-offset,coords_mm(:,2)-offset,coords_mm(:,3)];
            merstruct.currentmer.posterior.trajectory{side} = [trajectory(:,1)-offset,trajectory(:,2)-offset,trajectory(:,3)];
            merstruct.currentmer.medial.coords_mm{side} =  [coords_mm(:,1)-offset*2,coords_mm(:,2),coords_mm(:,3)]; %2mm medial (right is positive)
            merstruct.currentmer.medial.trajectory{side} = [trajectory(:,1)-offset*2,trajectory(:,2),trajectory(:,3)]; %2mm medial (right is positive)
        elseif side==2 %left
            merstruct.currentmer.central.coords_mm{side} = [coords_mm(:,1)+offset,coords_mm(:,2),coords_mm(:,3)];
            merstruct.currentmer.central.trajectory{side} = [trajectory(:,1)+offset,trajectory(:,2),trajectory(:,3)];
            merstruct.currentmer.anterior.coords_mm{side} = [coords_mm(:,1)+offset,coords_mm(:,2)+offset,coords_mm(:,3)];
            merstruct.currentmer.anterior.trajectory{side} = [trajectory(:,1)+offset,trajectory(:,2)+offset,trajectory(:,3)];
            merstruct.currentmer.posterior.coords_mm{side} = [coords_mm(:,1)+offset,coords_mm(:,2)-offset,coords_mm(:,3)];
            merstruct.currentmer.posterior.trajectory{side} = [trajectory(:,1)+offset,trajectory(:,2)-offset,trajectory(:,3)];
            merstruct.currentmer.medial.coords_mm{side} =  [coords_mm(:,1)+offset*2,coords_mm(:,2),coords_mm(:,3)]; %2mm medial (right is positive)
            merstruct.currentmer.medial.trajectory{side} = [trajectory(:,1)+offset*2,trajectory(:,2),trajectory(:,3)]; %2mm medial (right is positive)
        end
    case 'medial'
        merstruct.currentmer.medial.coords_mm{side} = coords_mm;
        merstruct.currentmer.medial.trajectory{side} = trajectory;
        if side==1 %right
            merstruct.currentmer.central.coords_mm{side} = [coords_mm(:,1)+offset,coords_mm(:,2),coords_mm(:,3)];
            merstruct.currentmer.central.trajectory{side} = [trajectory(:,1)+offset,trajectory(:,2),trajectory(:,3)];
            merstruct.currentmer.anterior.coords_mm{side} = [coords_mm(:,1)+offset,coords_mm(:,2)+offset,coords_mm(:,3)];
            merstruct.currentmer.anterior.trajectory{side} = [trajectory(:,1)+offset,trajectory(:,2)+offset,trajectory(:,3)];
            merstruct.currentmer.posterior.coords_mm{side} = [coords_mm(:,1)+offset,coords_mm(:,2)-offset,coords_mm(:,3)];
            merstruct.currentmer.posterior.trajectory{side} = [trajectory(:,1)+offset,trajectory(:,2)-offset,trajectory(:,3)];
            merstruct.currentmer.lateral.coords_mm{side} =  [coords_mm(:,1)+offset*2,coords_mm(:,2),coords_mm(:,3)]; %2mm medial (right is positive)
            merstruct.currentmer.lateral.trajectory{side} = [trajectory(:,1)+offset*2,trajectory(:,2),trajectory(:,3)]; %2mm medial (right is positive)
        elseif side==2 %left
            merstruct.currentmer.central.coords_mm{side} = [coords_mm(:,1)-offset,coords_mm(:,2),coords_mm(:,3)];
            merstruct.currentmer.central.trajectory{side} = [trajectory(:,1)-offset,trajectory(:,2),trajectory(:,3)];
            merstruct.currentmer.anterior.coords_mm{side} = [coords_mm(:,1)-offset,coords_mm(:,2)+offset,coords_mm(:,3)];
            merstruct.currentmer.anterior.trajectory{side} = [trajectory(:,1)-offset,trajectory(:,2)+offset,trajectory(:,3)];
            merstruct.currentmer.posterior.coords_mm{side} = [coords_mm(:,1)-offset,coords_mm(:,2)-offset,coords_mm(:,3)];
            merstruct.currentmer.posterior.trajectory{side} = [trajectory(:,1)-offset,trajectory(:,2)-offset,trajectory(:,3)];
            merstruct.currentmer.lateral.coords_mm{side} =  [coords_mm(:,1)-offset*2,coords_mm(:,2),coords_mm(:,3)]; %2mm medial (right is positive)
            merstruct.currentmer.lateral.trajectory{side} = [trajectory(:,1)-offset*2,trajectory(:,2),trajectory(:,3)]; %2mm medial (right is positive)
        end
end


function mer = ea_coords2mer(coords_mm,trajectory,resultfig,options)
        
offset = getfield(getappdata(resultfig,'merstruct'),'offset');
contact_length = getfield(getappdata(resultfig,'elspec'),'contact_length');

% x-axis --> negative = medial, positive = lateral
% y-axis --> negative = posterior, positive = anterior
% z-axis --> negative = inferior, positive = superior
for side=1:length(options.sides)
    dxyz = sqrt((diff(coords_mm{side}(1:2,1))^2)+(diff(coords_mm{side}(1:2,2))^2)+diff(coords_mm{side}(1:2,3))^2);
    slope = mean(diff(coords_mm{side}))/dxyz; %mean(diff(coords_mm{side}))/norm(mean(diff(coords_mm{side})))
    coords_mm{side} = coords_mm{side}-repmat(slope*contact_length/2,length(coords_mm{side}),1);
    trajectory{side} = trajectory{side}-repmat(slope*contact_length/2,length(trajectory{side}),1);
    if ea_getnativemni==1
        mer.central.coords_mm{side} = coords_mm{side};
        mer.central.trajectory{side} = trajectory{side};
        mer.anterior.coords_mm{side} = [coords_mm{side}(:,1),coords_mm{side}(:,2)+offset,coords_mm{side}(:,3)]; %2mm anterior
        mer.anterior.trajectory{side} = [trajectory{side}(:,1),trajectory{side}(:,2)+offset,trajectory{side}(:,3)]; %2mm anterior
        mer.posterior.coords_mm{side} = [coords_mm{side}(:,1),coords_mm{side}(:,2)-offset,coords_mm{side}(:,3)]; %2mm posterior
        mer.posterior.trajectory{side} = [trajectory{side}(:,1),trajectory{side}(:,2)-offset,trajectory{side}(:,3)]; %2mm posterior
        if side==1 %right
            mer.lateral.coords_mm{side} = [coords_mm{side}(:,1)+offset,coords_mm{side}(:,2),coords_mm{side}(:,3)]; %2mm lateral (right is positive)
            mer.lateral.trajectory{side} = [trajectory{side}(:,1)+offset,trajectory{side}(:,2),trajectory{side}(:,3)]; %2mm lateral (right is positive)
            mer.medial.coords_mm{side} = [coords_mm{side}(:,1)-offset,coords_mm{side}(:,2),coords_mm{side}(:,3)]; %2mm medial (right is positive)
            mer.medial.trajectory{side} = [trajectory{side}(:,1)-offset,trajectory{side}(:,2),trajectory{side}(:,3)]; %2mm medial (right is positive)
        elseif side==2 %left
            mer.lateral.coords_mm{side} = [coords_mm{side}(:,1)-offset,coords_mm{side}(:,2),coords_mm{side}(:,3)]; %2mm lateral (left is negative)
            mer.lateral.trajectory{side} = [trajectory{side}(:,1)-offset,trajectory{side}(:,2),trajectory{side}(:,3)]; %2mm lateral (left is negative)
            mer.medial.coords_mm{side} = [coords_mm{side}(:,1)+offset,coords_mm{side}(:,2),coords_mm{side}(:,3)]; %2mm medial (left is negative)
            mer.medial.trajectory{side} = [trajectory{side}(:,1)+offset,trajectory{side}(:,2),trajectory{side}(:,3)]; %2mm medial (left is negative)
        end
    end
    
    if ea_getnativemni==2
        mer.central.coords_mm{side} = coords_mm{side};
        mer.central.trajectory{side} = trajectory{side};
        mer.anterior.coords_mm{side} = [coords_mm{side}(:,1),coords_mm{side}(:,2)+offset,coords_mm{side}(:,3)]; %2mm anterior
        mer.anterior.trajectory{side} = [trajectory{side}(:,1),trajectory{side}(:,2)+offset,trajectory{side}(:,3)]; %2mm anterior
        mer.posterior.coords_mm{side} = [coords_mm{side}(:,1),coords_mm{side}(:,2)-offset,coords_mm{side}(:,3)]; %2mm posterior
        mer.posterior.trajectory{side} = [trajectory{side}(:,1),trajectory{side}(:,2)-offset,trajectory{side}(:,3)]; %2mm posterior
        if side==1 %right
            mer.lateral.coords_mm{side} = [coords_mm{side}(:,1)+offset,coords_mm{side}(:,2),coords_mm{side}(:,3)]; %2mm lateral (right is positive)
            mer.lateral.trajectory{side} = [trajectory{side}(:,1)+offset,trajectory{side}(:,2),trajectory{side}(:,3)]; %2mm lateral (right is positive)
            mer.medial.coords_mm{side} = [coords_mm{side}(:,1)-offset,coords_mm{side}(:,2),coords_mm{side}(:,3)]; %2mm medial (right is positive)
            mer.medial.trajectory{side} = [trajectory{side}(:,1)-offset,trajectory{side}(:,2),trajectory{side}(:,3)]; %2mm medial (right is positive)
        elseif side==2 %left
            mer.lateral.coords_mm{side} = [coords_mm{side}(:,1)-offset,coords_mm{side}(:,2),coords_mm{side}(:,3)]; %2mm lateral (left is negative)
            mer.lateral.trajectory{side} = [trajectory{side}(:,1)-offset,trajectory{side}(:,2),trajectory{side}(:,3)]; %2mm lateral (left is negative)
            mer.medial.coords_mm{side} = [coords_mm{side}(:,1)+offset,coords_mm{side}(:,2),coords_mm{side}(:,3)]; %2mm medial (left is negative)
            mer.medial.trajectory{side} = [trajectory{side}(:,1)+offset,trajectory{side}(:,2),trajectory{side}(:,3)]; %2mm medial (left is negative)
        end
    end
end
    


function refreshresultfig(handles,resultfig,options,side)
% this part makes changes of the figure active:
clearmertrajectories(handles,resultfig,options,side)
merstruct=getappdata(resultfig,'merstruct');
merhandles=getappdata(resultfig,'merhandles');
mermarkers=getappdata(resultfig,'mermarkers');
color=merstruct.colormap;
% Draw default trajectories
set(0,'CurrentFigure',resultfig); hold on;
for iSide = side
    
    if iSide==1 % right = 1 
        merhandles.central_right = plot3(merstruct.currentmer.central.trajectory{side}(:,1),merstruct.currentmer.central.trajectory{side}(:,2),merstruct.currentmer.central.trajectory{side}(:,3),'color',color(1,:),'linew',5,'tag','central_right');
        merhandles.anterior_right = plot3(merstruct.currentmer.anterior.trajectory{side}(:,1),merstruct.currentmer.anterior.trajectory{side}(:,2),merstruct.currentmer.anterior.trajectory{side}(:,3),'color',color(2,:),'linew',5,'tag','anterior_right');
        merhandles.posterior_right = plot3(merstruct.currentmer.posterior.trajectory{side}(:,1),merstruct.currentmer.posterior.trajectory{side}(:,2),merstruct.currentmer.posterior.trajectory{side}(:,3),'color',color(3,:),'linew',5,'tag','posterior_right');
        merhandles.lateral_right = plot3(merstruct.currentmer.lateral.trajectory{side}(:,1),merstruct.currentmer.lateral.trajectory{side}(:,2),merstruct.currentmer.lateral.trajectory{side}(:,3),'color',color(4,:),'linew',5,'tag','lateral_right');
        merhandles.medial_right = plot3(merstruct.currentmer.medial.trajectory{side}(:,1),merstruct.currentmer.medial.trajectory{side}(:,2),merstruct.currentmer.medial.trajectory{side}(:,3),'color',color(5,:),'linew',5,'tag','medial_right');
        set(handles.editimplanteddepth_right,'String','0')
        
    elseif iSide==2 % left = 2
        merhandles.central_left = plot3(merstruct.currentmer.central.trajectory{side}(:,1),merstruct.currentmer.central.trajectory{side}(:,2),merstruct.currentmer.central.trajectory{side}(:,3),'color',color(1,:),'linew',5,'tag','central_left');
        merhandles.anterior_left = plot3(merstruct.currentmer.anterior.trajectory{side}(:,1),merstruct.currentmer.anterior.trajectory{side}(:,2),merstruct.currentmer.anterior.trajectory{side}(:,3),'color',color(2,:),'linew',5,'tag','anterior_left');
        merhandles.posterior_left = plot3(merstruct.currentmer.posterior.trajectory{side}(:,1),merstruct.currentmer.posterior.trajectory{side}(:,2),merstruct.currentmer.posterior.trajectory{side}(:,3),'color',color(3,:),'linew',5,'tag','posterior_left');
        merhandles.lateral_left = plot3(merstruct.currentmer.lateral.trajectory{side}(:,1),merstruct.currentmer.lateral.trajectory{side}(:,2),merstruct.currentmer.lateral.trajectory{side}(:,3),'color',color(4,:),'linew',5,'tag','lateral_left');
        merhandles.medial_left = plot3(merstruct.currentmer.medial.trajectory{side}(:,1),merstruct.currentmer.medial.trajectory{side}(:,2),merstruct.currentmer.medial.trajectory{side}(:,3),'color',color(5,:),'linew',5,'tag','medial_left');    
        set(handles.editimplanteddepth_left,'String','0')
    end    
end

setappdata(resultfig,'merhandles',merhandles)


function clearkeyboardcontrolgroup(handles,mercontrolfig,resultfig,options)
set(handles.keycentral_left,'Value',0)
set(handles.keyanterior_left,'Value',0)
set(handles.keyposterior_left,'Value',0)
set(handles.keylateral_left,'Value',0)
set(handles.keymedial_left,'Value',0)

set(handles.keycentral_right,'Value',0)
set(handles.keyanterior_right,'Value',0)
set(handles.keyposterior_right,'Value',0)
set(handles.keycentral_right,'Value',0)
set(handles.keymedial_right,'Value',0)

setappdata(resultfig,'keymer',[])

function [files,subs] = fdir(Path,ext)
% Returns contents of path as struct
% Option to keep only files with input extension
%
% Inputs:
%
% 
% Ari Kappel, 2017
%

if isempty(Path) || ~exist('Path','var')
    Path = pwd;
end
    
contents = dir(Path);
contents = contents(cellfun(@(x) isempty(regexp(x, '^\.', 'once')), {contents.name}));

if ~isempty({contents(~[contents.isdir]).name})
    files = contents(~[contents.isdir]);
end
if ~isempty({contents([contents.isdir]).name})
    subs = contents([contents.isdir]);
end

if nargin==2
    files = files(~cellfun(@isempty , strfind({files.name},ext)));
end


function ea_keypress(handles,event)
% this is the main keypress function for the resultfigure. Add event
% listeners here.

if ismember('shift',event.Modifier)
end
