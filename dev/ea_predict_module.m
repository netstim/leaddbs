function varargout = ea_predict_module(varargin)
% EA_PREDICT_MODULE MATLAB code for ea_predict_module.fig
%      EA_PREDICT_MODULE, by itself, creates a new EA_PREDICT_MODULE or raises the existing
%      singleton*.
%
%      H = EA_PREDICT_MODULE returns the handle to a new EA_PREDICT_MODULE or the handle to
%      the existing singleton*.
%
%      EA_PREDICT_MODULE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in EA_PREDICT_MODULE.M with the given input arguments.
%
%      EA_PREDICT_MODULE('Property','Value',...) creates a new EA_PREDICT_MODULE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ea_predict_module_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ea_predict_module_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ea_predict_module

% Last Modified by GUIDE v2.5 20-Feb-2017 13:43:13

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ea_predict_module_OpeningFcn, ...
                   'gui_OutputFcn',  @ea_predict_module_OutputFcn, ...
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


% --- Executes just before ea_predict_module is made visible.
function ea_predict_module_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ea_predict_module (see VARARGIN)

% Choose default command line output for ea_predict_module
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ea_predict_module wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = ea_predict_module_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
