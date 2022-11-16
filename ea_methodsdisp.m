function varargout = ea_methodsdisp(varargin)
% EA_METHODSDISP MATLAB code for ea_methodsdisp.fig
% This will show the method window.
% Requires a c formatted string (printf) as input, within a cell array
% String passed requires to be terminated with '\n'
% Text before the 5th character will be ignored
% Example: {'\n\nText to print\nReferences: 1) ref1\n'}
% will appear as :
%   Text to print
%   References: 1) ref1
%
%
%      EA_METHODSDISP, by itself, creates a new EA_METHODSDISP or raises the existing
%      singleton*.
%
%      H = EA_METHODSDISP returns the handle to a new EA_METHODSDISP or the handle to
%      the existing singleton*.
%
%      EA_METHODSDISP('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in EA_METHODSDISP.M with the given input arguments.
%
%      EA_METHODSDISP('Property','Value',...) creates a new EA_METHODSDISP or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ea_methodsdisp_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ea_methodsdisp_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ea_methodsdisp

% Last Modified by GUIDE v2.5 19-Mar-2017 16:58:04

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ea_methodsdisp_OpeningFcn, ...
                   'gui_OutputFcn',  @ea_methodsdisp_OutputFcn, ...
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


% --- Executes just before ea_methodsdisp is made visible.
function ea_methodsdisp_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ea_methodsdisp (see VARARGIN)


if isempty(handles.methodstxt.String)
    %this will parse the string that is passed as c (printf) formatted string
    handles.methodstxt.String=ea_strsplit(sprintf(varargin{1}{1}(5:end)));
else
    try
        %this will parse the string that is passed as c (printf) formatted string
        handles.methodstxt.String=ea_strsplit([sprintf('%s\n\n',handles.methodstxt.String),sprintf(varargin{1}{1})]);
    catch
        %this will parse the string that is passed as c (printf) formatted string
        handles.methodstxt.String=ea_strsplit(sprintf(varargin{1}{1}(5:end)));
    end
end

set(hObject,'Name','Used Methods');

movegui(hObject,'northeast');
% Choose default command line output for ea_methodsdisp
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ea_methodsdisp wait for user response (see UIRESUME)
% uiwait(handles.methodsfig);

function cellstr=ea_strsplit(str)
cnt=1;
ns=strfind(str,sprintf('\n'));
offs=1;

for n=1:length(ns)
    cellstr{cnt}=str(offs:ns(n)-1);
    offs=ns(n)+1;
    cnt=cnt+1;
end


% --- Outputs from this function are returned to the command line.
function varargout = ea_methodsdisp_OutputFcn(hObject, eventdata, handles) 
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