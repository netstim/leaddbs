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

% Last Modified by GUIDE v2.5 15-Jan-2020 17:13:07

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
movegui(hObject,'northwest');



setappdata(handles.trajectorycontrol,'chandles',handles);
obj=varargin{1};
% add menu:
f = uimenu('Label','Tools');
uimenu(f,'Label','Export Plan as Reconstruction...','Callback',{@ea_plan2reconstruction,obj},'Accelerator','E');

set(handles.electrode_model_plan,'String',ea_resolve_elspec);

setappdata(obj.plotFigureH,'trajcontrolfig',handles.trajectorycontrol);
setappdata(handles.trajectorycontrol,'obj',obj);
set(handles.trajectorycontrol,'name','Edit Trajectory');
ea_synctrajectoryhandles(handles,obj)


function ea_plan2reconstruction(~,~,obj)
[FileName,PathName] = uiputfile('ea_reconstruction.mat','Choose destination of ea_reconstruction.mat');

if ~strcmp(FileName,'ea_reconstruction.mat')
    ea_warning('Files that are named differently than ea_reconstruction.mat will not be recognized by Lead-DBS.');
end
elstruct=obj.plan2elstruct;
options=obj.options;
[options.root,options.patientname]=fileparts(PathName);
ea_save_reconstruction(elstruct(1).coords_mm,elstruct(1).trajectory,elstruct(1).markers,obj.plan2elstruct_model,1,options,FileName);

function ea_synctrajectoryhandles(handles,obj)
% set handles to match obj

% set coordinates
if obj.hasPlanning
set(handles.targetX,'String',num2str(obj.target.target(1))); set(handles.targetY,'String',num2str(obj.target.target(2))); set(handles.targetZ,'String',num2str(obj.target.target(3)));
set(handles.entryX,'String',num2str(obj.target.entry(1))); set(handles.entryY,'String',num2str(obj.target.entry(2))); set(handles.entryZ,'String',num2str(obj.target.entry(3)));
end
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

%% set showplanning.
set(handles.showPlanning,'Value',obj.showPlanning);
set(handles.showPlanning,'enable',ea_bool2onoff(obj.hasPlanning));
% subordinate enables
if ~get(handles.showPlanning,'Value') || ~ea_bool2onoff(get(handles.showPlanning,'enable'))
    onoff='off';
else
    onoff='on';
end
set(handles.targetX,'enable',onoff); set(handles.targetY,'enable',onoff); set(handles.targetZ,'enable',onoff);
set(handles.entryX,'enable',onoff); set(handles.entryY,'enable',onoff); set(handles.entryZ,'enable',onoff);
set(handles.space,'enable',onoff);
set(handles.color,'enable',onoff);
switch obj.planningAppearance
    case 'line'
        set(handles.planningappearance,'Value',1);
        set(handles.electrode_model_plan,'Enable','off');
    case 'electrode'
        set(handles.planningappearance,'Value',2);
        set(handles.electrode_model_plan,'Enable','on');
end

string_list = get(handles.electrode_model_plan,'String');
[~,whichentry]=ismember(obj.plan2elstruct_model,string_list);
set(handles.electrode_model_plan,'Value',whichentry);


if ~(get(handles.showPlanning,'Value')) || ~ea_bool2onoff(get(handles.showPlanning,'enable')) || get(handles.space,'Value')>1
    onoff='off';
else
    onoff='on';
end
set(handles.right,'enable',onoff);
set(handles.left,'enable',onoff);
set(handles.anterior,'enable',onoff);
set(handles.posterior,'enable',onoff);
set(handles.ventral,'enable',onoff);
set(handles.dorsal,'enable',onoff);
set(handles.AC,'enable',onoff);
set(handles.PC,'enable',onoff);
set(handles.MCP,'enable',onoff);

%% macro/DBS electrode
set(handles.showMacro,'Value',obj.showMacro);
set(handles.showMacro,'enable',ea_bool2onoff(obj.hasMacro));
set(handles.electrode_model_popup,'String',ea_resolve_elspec);
% subordinate enables
if ~(get(handles.showMacro,'Value')) || ~ea_bool2onoff(get(handles.showMacro,'enable'))
    onoff='off';
else
    onoff='on';
end
set(handles.electrode_model_popup,'enable',onoff);

%% micro/MER
set(handles.showMicro,'Value',obj.showMicro);
set(handles.relateMicro,'String',{'Base location on Macroelectrode','Base location on Microelectrode'});
switch obj.relateMicro
    case 'macro'
        set(handles.relateMicro,'Value',1);
    case 'planning'
        set(handles.relateMicro,'Value',2);
end
% subordinate enables
set(handles.relateMicro,'enable',ea_bool2onoff(get(handles.showMicro,'Value')));

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
obj = getappdata(handles.trajectorycontrol,'obj');
obj.showPlanning = get(handles.showPlanning,'Value');
obj.togglestates(1) = get(handles.showPlanning,'Value');
obj.toggleH.State = get(handles.showPlanning,'Value');
ea_synctrajectoryhandles(handles,obj);

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
c = ea_uisetcolor;
if any(c)
obj=getappdata(handles.trajectorycontrol,'obj');
obj.color=c;
ea_synctrajectoryhandles(handles,obj);
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
else
    set(hObject,'Value',1);
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
else
    set(hObject,'Value',1);
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
else
    set(hObject,'Value',1);
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
ea_synctrajectoryhandles(handles,obj);

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
else
    set(hObject,'Value',1);
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
    else
    set(hObject,'Value',1);
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
    else
    set(hObject,'Value',1);
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
    else
    set(hObject,'Value',1);
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
    else
    set(hObject,'Value',1);
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
    else
    set(hObject,'Value',1);
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


% --- Executes on button press in showMacro.
function showMacro_Callback(hObject, eventdata, handles)
% hObject    handle to showMacro (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of showMacro
obj=getappdata(handles.trajectorycontrol,'obj');
obj.showMacro=get(hObject,'Value');
obj.togglestates(2)=get(handles.showMacro,'Value');
ea_synctrajectoryhandles(handles,obj);

% --- Executes on selection change in plan_electrode_model.
function plan_electrode_model_Callback(hObject, eventdata, handles)
% hObject    handle to plan_electrode_model (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns plan_electrode_model contents as cell array
%        contents{get(hObject,'Value')} returns selected item from plan_electrode_model
obj=getappdata(handles.trajectorycontrol,'obj');
    options.elmodeln = get(handles.plan_electrode_model,'Value');
    string_list = get(handles.plan_electrode_model,'String');
    obj.elmodel=string_list{options.elmodeln};
ea_synctrajectoryhandles(handles,obj);


% --- Executes during object creation, after setting all properties.
function plan_electrode_model_CreateFcn(hObject, eventdata, handles)
% hObject    handle to plan_electrode_model (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in showMicro.
function showMicro_Callback(hObject, eventdata, handles)
% hObject    handle to showMicro (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of showMicro
obj=getappdata(handles.trajectorycontrol,'obj');
obj.showMicro=get(hObject,'Value');
obj.togglestates(3)=get(handles.showMicro,'Value');
ea_synctrajectoryhandles(handles,obj);


% --- Executes on selection change in relateMicro.
function relateMicro_Callback(hObject, eventdata, handles)
% hObject    handle to relateMicro (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns relateMicro contents as cell array
%        contents{get(hObject,'Value')} returns selected item from relateMicro
obj=getappdata(handles.trajectorycontrol,'obj');
switch get(hObject,'Value')
    case 1
        obj.relateMicro='macro';
    case 2
        obj.relateMicro='planning';
end
ea_synctrajectoryhandles(handles,obj);


% --- Executes during object creation, after setting all properties.
function relateMicro_CreateFcn(hObject, eventdata, handles)
% hObject    handle to relateMicro (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in addtraj.
function addtraj_Callback(hObject, eventdata, handles)
% hObject    handle to addtraj (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
obj=getappdata(handles.trajectorycontrol,'obj');
if obj.hasPlanning % -> This will add a new trajectory unrelated to the present one


else

    obj.target=ea_getstandardtarget(obj.side);
    obj.showPlanning=1;
    obj.hasPlanning=1;
    obj.showPlanning=1;
    ea_synctrajectoryhandles(handles,obj);

end


% --- Executes when user attempts to close trajectorycontrol.
function trajectorycontrol_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to trajectorycontrol (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
obj=getappdata(handles.trajectorycontrol,'obj');
delete(hObject);

ea_save_trajectory(obj);


% --- Executes on selection change in planningappearance.
function planningappearance_Callback(hObject, eventdata, handles)
% hObject    handle to planningappearance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns planningappearance contents as cell array
%        contents{get(hObject,'Value')} returns selected item from planningappearance
obj=getappdata(handles.trajectorycontrol,'obj');
switch get(hObject,'value')
    case 1
        obj.planningAppearance='line';
    case 2
        obj.planningAppearance='electrode';
end
ea_synctrajectoryhandles(handles,obj);


% --- Executes during object creation, after setting all properties.
function planningappearance_CreateFcn(hObject, eventdata, handles)
% hObject    handle to planningappearance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in electrode_model_popup.
function electrode_model_popup_Callback(hObject, eventdata, handles)
% hObject    handle to electrode_model_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns electrode_model_popup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from electrode_model_popup
obj=getappdata(handles.trajectorycontrol,'obj');
    options.elmodeln = get(handles.electrode_model_popup,'Value');
    string_list = get(handles.electrode_model_popup,'String');
    obj.elmodel=string_list{options.elmodeln};
ea_synctrajectoryhandles(handles,obj);

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


% --- Executes on selection change in electrode_model_plan.
function electrode_model_plan_Callback(hObject, eventdata, handles)
% hObject    handle to electrode_model_plan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns electrode_model_plan contents as cell array
%        contents{get(hObject,'Value')} returns selected item from electrode_model_plan
obj=getappdata(handles.trajectorycontrol,'obj');
    options.elmodeln = get(handles.electrode_model_plan,'Value');
    string_list = get(handles.electrode_model_plan,'String');
    obj.plan2elstruct_model=string_list{options.elmodeln};
ea_synctrajectoryhandles(handles,obj);

% --- Executes during object creation, after setting all properties.
function electrode_model_plan_CreateFcn(hObject, eventdata, handles)
% hObject    handle to electrode_model_plan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
