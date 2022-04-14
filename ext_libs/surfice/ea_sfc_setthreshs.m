function varargout = ea_sfc_setthreshs(varargin)
% EA_SFC_SETTHRESHS MATLAB code for ea_sfc_setthreshs.fig
%      EA_SFC_SETTHRESHS, by itself, creates a new EA_SFC_SETTHRESHS or raises the existing
%      singleton*.
%
%      H = EA_SFC_SETTHRESHS returns the handle to a new EA_SFC_SETTHRESHS or the handle to
%      the existing singleton*.
%
%      EA_SFC_SETTHRESHS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in EA_SFC_SETTHRESHS.M with the given input arguments.
%
%      EA_SFC_SETTHRESHS('Property','Value',...) creates a new EA_SFC_SETTHRESHS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ea_sfc_setthreshs_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ea_sfc_setthreshs_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ea_sfc_setthreshs

% Last Modified by GUIDE v2.5 17-Mar-2017 12:53:27

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ea_sfc_setthreshs_OpeningFcn, ...
                   'gui_OutputFcn',  @ea_sfc_setthreshs_OutputFcn, ...
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


% --- Executes just before ea_sfc_setthreshs is made visible.
function ea_sfc_setthreshs_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ea_sfc_setthreshs (see VARARGIN)

% Choose default command line output for ea_sfc_setthreshs
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);


threshs=varargin{1};
fnames=varargin{2};
setappdata(handles.threshsfigure,'threshs',threshs);
handles.threshs.Data=threshs;
handles.threshs.ColumnEditable=logical([1 1 1 1]);
handles.threshs.ColumnName={'MinPos','MaxPos','MaxNeg','MinNeg'};
handles.threshs.RowName=fnames;

% UIWAIT makes ea_sfc_setthreshs wait for user response (see UIRESUME)
 uiwait(handles.threshsfigure);


% --- Outputs from this function are returned to the command line.
function varargout = ea_sfc_setthreshs_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure

varargout{1} = getappdata(handles.threshsfigure,'action');
varargout{2} = getappdata(handles.threshsfigure,'threshs');

% The figure can be deleted now
delete(hObject);

% --- Executes on button press in proceedbutton.
function proceedbutton_Callback(hObject, eventdata, handles)
% hObject    handle to proceedbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
setappdata(handles.threshsfigure,'action','proceed');
setappdata(handles.threshsfigure,'threshs',handles.threshs.Data);

threshsfigure_CloseRequestFcn(handles.threshsfigure,[],handles);

% --- Executes on button press in cancelbutton.
function cancelbutton_Callback(hObject, eventdata, handles)
% hObject    handle to cancelbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
setappdata(handles.threshsfigure,'action','cancel');
threshsfigure_CloseRequestFcn(handles.threshsfigure,[],handles);

% --- Executes on button press in propagatebutton.
function propagatebutton_Callback(hObject, eventdata, handles)
% hObject    handle to propagatebutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
threshs=handles.threshs.Data;
threshs=repmat(threshs(1,:),size(threshs,1),1);
handles.threshs.Data=threshs;

% --- Executes when user attempts to close threshsfigure.
function threshsfigure_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to threshsfigure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
%delete(hObject);

if isequal(get(hObject, 'waitstatus'), 'waiting')
% The GUI is still in UIWAIT, us UIRESUME
uiresume(hObject);
else
% The GUI is no longer waiting, just close it
delete(hObject);
end
