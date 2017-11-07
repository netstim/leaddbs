function varargout = ea_trajectorycontrol(varargin)
% EA_TRAJECTORYCONTROL MATLAB code for ea_trajectorycontrol.fig
%      EA_TRAJECTORYCONTROL, by itself, creates a new EA_TRAJECTORYCONTROL or raises the existing
%      singleton*.
%
%      H = EA_TRAJECTORYCONTROL returns the handle to a new EA_TRAJECTORYCONTROL or the handle to
%      the existing singleton*.
%
%      EA_TRAJECTORYCONTROL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in EA_TRAJECTORYCONTROL.M with the given input arguments.
%
%      EA_TRAJECTORYCONTROL('Property','Value',...) creates a new EA_TRAJECTORYCONTROL or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ea_trajectorycontrol_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ea_trajectorycontrol_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ea_trajectorycontrol

% Last Modified by GUIDE v2.5 05-Nov-2017 13:39:06

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ea_trajectorycontrol_OpeningFcn, ...
                   'gui_OutputFcn',  @ea_trajectorycontrol_OutputFcn, ...
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


% --- Executes just before ea_trajectorycontrol is made visible.
function ea_trajectorycontrol_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ea_trajectorycontrol (see VARARGIN)

% Choose default command line output for ea_trajectorycontrol
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ea_trajectorycontrol wait for user response (see UIRESUME)
% uiwait(handles.trajectorycontrol);

setappdata(handles.trajectorycontrol,'chandles',handles);
obj=varargin{1};
setappdata(handles.trajectorycontrol,'obj',obj);
ea_synctrajectoryhandles(handles,obj)

function ea_synctrajectoryhandles(handles,obj)
% set handles to match obj

% set coordinates
set(handles.targetX,'String',num2str(obj.target.target(1))); set(handles.targetY,'String',num2str(obj.target.target(2))); set(handles.targetZ,'String',num2str(obj.target.target(3)));
set(handles.entryX,'String',num2str(obj.target.entry(1))); set(handles.entryY,'String',num2str(obj.target.entry(2))); set(handles.entryZ,'String',num2str(obj.target.entry(3)));

% set space
set(handles.space,'Value',obj.planRelative(5));

% set relative to radio buttons
switch obj.planRelative(1)
    case 1
        set(handles.AC,'Value',1);
    case 2
        set(handles.MCP,'Value',1);
    case 3
        set(handles.PC,'Value',1);
end
switch obj.planRelative(2)
    case 1
        set(handles.right,'Value',1);
    case 2
        set(handles.left,'Value',1);
end
switch obj.planRelative(3)
    case 1
        set(handles.anterior,'Value',1);
    case 2
        set(handles.posterior,'Value',1);
end
switch obj.planRelative(4)
    case 1
        set(handles.ventral,'Value',1);
    case 2
        set(handles.dorsal,'Value',1);
end

% set color backgroundcolor
set(handles.color,'BackgroundColor',obj.color);
% set showplanning.
set(handles.showPlanning,'Value',obj.showPlanning);

% --- Outputs from this function are returned to the command line.
function varargout = ea_trajectorycontrol_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in showPlanning.
function showPlanning_Callback(hObject, eventdata, handles)
% hObject    handle to showPlanning (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of showPlanning
obj=getappdata(handles.trajectorycontrol,'obj');
obj.showPlanning=get(handles.showPlanning,'Value');


function targetX_Callback(hObject, eventdata, handles)
% hObject    handle to targetX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of targetX as text
%        str2double(get(hObject,'String')) returns contents of targetX as a double
obj=getappdata(handles.trajectorycontrol,'obj');
obj.target.target(1)=str2double(get(handles.targetX,'String'));


% --- Executes during object creation, after setting all properties.
function targetX_CreateFcn(hObject, eventdata, handles)
% hObject    handle to targetX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function targetY_Callback(hObject, eventdata, handles)
% hObject    handle to txtY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtY as text
%        str2double(get(hObject,'String')) returns contents of txtY as a double
obj=getappdata(handles.trajectorycontrol,'obj');
obj.target.target(2)=str2double(get(handles.targetY,'String'));

% --- Executes during object creation, after setting all properties.
function txtY_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function targetZ_Callback(hObject, eventdata, handles)
% hObject    handle to targetZ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of targetZ as text
%        str2double(get(hObject,'String')) returns contents of targetZ as a double
obj=getappdata(handles.trajectorycontrol,'obj');
obj.target.target(3)=str2double(get(handles.targetZ,'String'));

% --- Executes during object creation, after setting all properties.
function targetZ_CreateFcn(hObject, eventdata, handles)
% hObject    handle to targetZ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function entryX_Callback(hObject, eventdata, handles)
% hObject    handle to entryX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of entryX as text
%        str2double(get(hObject,'String')) returns contents of entryX as a double
obj=getappdata(handles.trajectorycontrol,'obj');
obj.target.entry(1)=str2double(get(handles.entryX,'String'));

% --- Executes during object creation, after setting all properties.
function entryX_CreateFcn(hObject, eventdata, handles)
% hObject    handle to entryX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function entryY_Callback(hObject, eventdata, handles)
% hObject    handle to entryY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of entryY as text
%        str2double(get(hObject,'String')) returns contents of entryY as a double
obj=getappdata(handles.trajectorycontrol,'obj');
obj.target.entry(2)=str2double(get(handles.entryY,'String'));

% --- Executes during object creation, after setting all properties.
function entryY_CreateFcn(hObject, eventdata, handles)
% hObject    handle to entryY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function entryZ_Callback(hObject, eventdata, handles)
% hObject    handle to entryZ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of entryZ as text
%        str2double(get(hObject,'String')) returns contents of entryZ as a double
obj=getappdata(handles.trajectorycontrol,'obj');
obj.target.entry(3)=str2double(get(handles.entryZ,'String'));

% --- Executes during object creation, after setting all properties.
function entryZ_CreateFcn(hObject, eventdata, handles)
% hObject    handle to entryZ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in color.
function color_Callback(hObject, eventdata, handles)
% hObject    handle to color (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
c=uisetcolor;
if any(c)
obj=getappdata(handles.trajectorycontrol,'obj');
obj.color=c;
set(handles.color,'BackgroundColor',c);
end

% --- Executes on button press in AC.
function AC_Callback(hObject, eventdata, handles)
% hObject    handle to AC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of AC
if get(hObject,'Value')
    set(handles.MCP,'Value',0);
    set(handles.PC,'Value',0);
    obj=getappdata(handles.trajectorycontrol,'obj');
    obj.planRelative(1)=1;
end

% --- Executes on button press in MCP.
function MCP_Callback(hObject, eventdata, handles)
% hObject    handle to MCP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of MCP
if get(hObject,'Value')
    set(handles.AC,'Value',0);
    set(handles.PC,'Value',0);
    obj=getappdata(handles.trajectorycontrol,'obj');
    obj.planRelative(1)=2;
end

% --- Executes on button press in PC.
function PC_Callback(hObject, eventdata, handles)
% hObject    handle to PC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of PC
if get(hObject,'Value')
    set(handles.MCP,'Value',0);
    set(handles.AC,'Value',0);
    obj=getappdata(handles.trajectorycontrol,'obj');
    obj.planRelative(1)=3;
end

% --- Executes on selection change in space.
function space_Callback(hObject, eventdata, handles)
% hObject    handle to space (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns space contents as cell array
%        contents{get(hObject,'Value')} returns selected item from space

    obj=getappdata(handles.trajectorycontrol,'obj');
    obj.planRelative(5)=get(handles.space,'Value');

% --- Executes during object creation, after setting all properties.
function space_CreateFcn(hObject, eventdata, handles)
% hObject    handle to space (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in right.
function right_Callback(hObject, eventdata, handles)
% hObject    handle to right (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of right
if get(hObject,'Value')
    set(handles.left,'Value',0);
    obj=getappdata(handles.trajectorycontrol,'obj');
    obj.planRelative(2)=1;
end

% --- Executes on button press in left.
function left_Callback(hObject, eventdata, handles)
% hObject    handle to left (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of left
if get(hObject,'Value')
    set(handles.right,'Value',0);
    obj=getappdata(handles.trajectorycontrol,'obj');
    obj.planRelative(2)=2;
end

% --- Executes on button press in anterior.
function anterior_Callback(hObject, eventdata, handles)
% hObject    handle to anterior (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of anterior
if get(hObject,'Value')
    set(handles.posterior,'Value',0);
    obj=getappdata(handles.trajectorycontrol,'obj');
    obj.planRelative(3)=1;
end

% --- Executes on button press in posterior.
function posterior_Callback(hObject, eventdata, handles)
% hObject    handle to posterior (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of posterior
if get(hObject,'Value')
    set(handles.anterior,'Value',0);
    obj=getappdata(handles.trajectorycontrol,'obj');
    obj.planRelative(3)=2;
end

% --- Executes on button press in ventral.
function ventral_Callback(hObject, eventdata, handles)
% hObject    handle to ventral (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ventral
if get(hObject,'Value')
    set(handles.dorsal,'Value',0);
    obj=getappdata(handles.trajectorycontrol,'obj');
    obj.planRelative(4)=1;
end

% --- Executes on button press in dorsal.
function dorsal_Callback(hObject, eventdata, handles)
% hObject    handle to dorsal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of dorsal
if get(hObject,'Value')
    set(handles.ventral,'Value',0);
    obj=getappdata(handles.trajectorycontrol,'obj');
    obj.planRelative(4)=2;
end


% --- Executes during object creation, after setting all properties.
function targetY_CreateFcn(hObject, eventdata, handles)
% hObject    handle to targetY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
