function varargout = ea_textdisp(varargin)
% EA_TEXTDISP MATLAB code for ea_textdisp.fig
%      EA_TEXTDISP, by itself, creates a new EA_TEXTDISP or raises the existing
%      singleton*.
%
%      H = EA_TEXTDISP returns the handle to a new EA_TEXTDISP or the handle to
%      the existing singleton*.
%
%      EA_TEXTDISP('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in EA_TEXTDISP.M with the given input arguments.
%
%      EA_TEXTDISP('Property','Value',...) creates a new EA_TEXTDISP or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ea_textdisp_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ea_textdisp_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ea_textdisp

% Last Modified by GUIDE v2.5 04-Apr-2017 11:16:27

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ea_textdisp_OpeningFcn, ...
                   'gui_OutputFcn',  @ea_textdisp_OutputFcn, ...
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


% --- Executes just before ea_textdisp is made visible.
function ea_textdisp_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ea_textdisp (see VARARGIN)


    handles.methodstxt.String=ea_strsplit(sprintf('%s',varargin{1}{1}));


set(hObject,'Name','Processing Report');

movegui(hObject,'northwest');
% Choose default command line output for ea_textdisp
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ea_textdisp wait for user response (see UIRESUME)
% uiwait(handles.methodsfig);

function cellstr=ea_strsplit(str)
cnt=1;
ns=strfind(str,'\n');
offs=1;

for n=1:length(ns)
    cellstr{cnt}=str(offs:ns(n)-1);
    offs=ns(n)+2;
    cnt=cnt+1;
end


% --- Outputs from this function are returned to the command line.
function varargout = ea_textdisp_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function methodstxt_Callback(hObject, eventdata, handles)
% hObject    handle to methodstxt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of methodstxt as text
%        str2double(get(hObject,'String')) returns contents of methodstxt as a double


% --- Executes during object creation, after setting all properties.
function methodstxt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to methodstxt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
