function varargout = editSessionGUI(varargin)
% EDITSESSIONGUI M-file for editSessionGUI.fig
%      EDITSESSIONGUI, by itself, creates a new EDITSESSIONGUI or raises the existing
%      singleton*.
%
%      H = EDITSESSIONGUI returns the handle to a new EDITSESSIONGUI or the handle to
%      the existing singleton*.
%
%      EDITSESSIONGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in EDITSESSIONGUI.M with the given input arguments.
%
%      EDITSESSIONGUI('Property','Value',...) creates a new EDITSESSIONGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before editSessionGUI_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to editSessionGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Copyright 2002-2003 The MathWorks, Inc.

% Edit the above text to modify the response to help editSessionGUI

% Last Modified by GUIDE v2.5 15-May-2005 12:33:07

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @editSessionGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @editSessionGUI_OutputFcn, ...
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


% --- Executes just before editSessionGUI is made visible.
function editSessionGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to editSessionGUI (see VARARGIN)

% Initialize session
% editSessionGUI can be called by passing an initial session
% structure:
%     session = editSessionGUI('session',session);
% or by passing nothing in which case the structure is built from scratch
%     session = editSessionGUI;
if length(varargin) == 2
    field = varargin{1};
    val = varargin{2};
    switch field
        case 'session'
            session = val;
        otherwise
            warn('editSessionGUI: invalid initialization argument')
    end
end

if ~exist('session','var')
    session = [];
end
if ~isfield(session,'mrLoadRetVersion')
    session.mrLoadRetVersion = mrLoadRetVersion;
end
if ~isfield(session,'description')
    session.description = '';
end
if ~isfield(session,'subject')
    session.subject = '';
end
if ~isfield(session,'operator')
    session.operator = '';
end
if ~isfield(session,'magnet')
    str = get(handles.magnetPopup,'string');
    session.magnet = str{1};
end
if ~isfield(session,'coil')
    str = get(handles.coilPopup,'string');
    session.coil = str{1};
end
if ~isfield(session,'protocol')
    session.protocol = '';
end

% Initialize gui from session
set(handles.descriptionText,'string',session.description);
set(handles.subjectText,'string',session.subject);
set(handles.operatorText,'string',session.operator);
set(handles.protocolText,'string',session.protocol);
magnetValue = find(strcmp(session.magnet,...
    get(handles.magnetPopup,'String')));
set(handles.magnetPopup,'Value',magnetValue);
coilValue = find(strcmp(session.coil,...
    get(handles.coilPopup,'String')));
set(handles.coilPopup,'Value',coilValue);

% Update handles structure
handles.session = session;
guidata(hObject, handles);

% Choose default command line output for editSessionGUI
handles.output = handles.session;

% UIWAIT makes editSessionGUI wait for user response (see UIRESUME)
uiwait(handles.figure);


% --- Outputs from this function are returned to the command line.
function varargout = editSessionGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
delete(handles.figure);


% --- descriptionText
function descriptionText_Callback(hObject, eventdata, handles)
% Update session.description
handles.session.description = get(hObject,'string');
guidata(hObject, handles);

function descriptionText_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- subjectText
function subjectText_Callback(hObject, eventdata, handles)
% Update session.subject
handles.session.subject = get(hObject,'string');
guidata(hObject, handles);

function subjectText_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- operatorText
function operatorText_Callback(hObject, eventdata, handles)
% Update session.operator
handles.session.operator = get(hObject,'string');
guidata(hObject, handles);

function operatorText_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- protocolText
function protocolText_Callback(hObject, eventdata, handles)
% Update session.protocol
handles.session.protocol = get(hObject,'string');
guidata(hObject, handles);

function protocolText_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- magnetPopup
function magnetPopup_Callback(hObject, eventdata, handles)
% Update session.magnet
str = get(hObject,'string');
val = get(hObject,'value');
handles.session.magnet = str{val};
guidata(hObject, handles);

function magnetPopup_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- coilPopup
function coilPopup_Callback(hObject, eventdata, handles)
% Update session.coil
str = get(hObject,'string');
val = get(hObject,'value');
handles.session.coil = str{val};
guidata(hObject, handles);

function coilPopup_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- okButton
function okButton_Callback(hObject, eventdata, handles)
handles.output = handles.session;
guidata(hObject, handles);
uiresume;


% --- cancelButton
function cancelButton_Callback(hObject, eventdata, handles)
handles.output = [];
guidata(hObject, handles);
uiresume;



