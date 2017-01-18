function varargout = lead(varargin)
% LEAD MATLAB code for lead.fig
%      LEAD, by itself, creates a new LEAD or raises the existing
%      singleton*.
%
%      H = LEAD returns the handle to a new LEAD or the handle to
%      the existing singleton*.
%
%      LEAD('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in LEAD.M with the given input arguments.
%
%      LEAD('Property','Value',...) creates a new LEAD or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before lead_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to lead_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help lead

% Last Modified by GUIDE v2.5 31-Jul-2016 21:00:09

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @lead_OpeningFcn, ...
    'gui_OutputFcn',  @lead_OutputFcn, ...
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


% --- Executes just before lead is made visible.
function lead_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to lead (see VARARGIN)

% Choose default command line output for lead
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
 try
    ea_compat_data; 
 end
earoot=ea_getearoot;

ea_setpath;
ea_checkbuildspace;

ea_checkspm;

% check existence of directories

ea_checkleaddirs;


% check for commands first
if nargin>3
    switch varargin{1}
        case 'dbs'
            lead_dbs;
            delete(handles.leadfigure)
            return
        case 'group';
            lead_group;
            delete(handles.leadfigure)
            return
        case 'connectome';
            lead_connectome;
            delete(handles.leadfigure)
            return
        case {'mapper'};
            lead_mapper;
            delete(handles.leadfigure)
            return
        case 'anatomy';
            lead_anatomy;
            delete(handles.leadfigure)
            return
        case 'macaque';
            lead_dbs macaque;
            delete(handles.leadfigure)
            return
        case 'version'
            disp(ea_getvsn('local'));
            delete(handles.leadfigure)
            return
        case 'speak'
            fprintf('\n \n \n \n %s \n \n','L337-D8Z: "H3LLo 7o joO MY Phr13nd. l1V3 Lon9 4nd pRO5P3r."'); % yes, this indeed is an easter-egg.
            delete(handles.leadfigure)
            return
    end
    
end


set(handles.leadfigure,'name','Welcome to the Lead Neuroimaging Suite');


set(0,'CurrentFigure',handles.leadfigure);
im=imread([earoot,'icons',filesep,'logo_lead.png']);
image(im);
axis off;
axis equal;

% add logos

ea_setbuttonbackdrop(handles.startdbs,[earoot,'icons',filesep,'logo_lead_dbs_small.png']);
ea_setbuttonbackdrop(handles.startconnectome,[earoot,'icons',filesep,'logo_lead_connectome_small.png']);
ea_setbuttonbackdrop(handles.startgroup,[earoot,'icons',filesep,'logo_lead_group_small.png']);
ea_setbuttonbackdrop(handles.startanatomy,[earoot,'icons',filesep,'logo_lead_anatomy_small.png']);
set(handles.versiontxt,'String',['v',ea_getvsn('local')]);




% UIWAIT makes lead wait for user response (see UIRESUME)
% uiwait(handles.leadfigure);


function ea_setbuttonbackdrop(buttonhandle,im)
im=imread(im);

set(buttonhandle,'CData',im);
set(buttonhandle,'String','');


% --- Outputs from this function are returned to the command line.
function varargout = lead_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure


% --- Executes on button press in startdbs.
function startdbs_Callback(hObject, eventdata, handles)
% hObject    handle to startdbs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

lead dbs;

% --- Executes on button press in startconnectome.
function startconnectome_Callback(hObject, eventdata, handles)
% hObject    handle to startconnectome (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
lead connectome;

% --- Executes on button press in startgroup.
function startgroup_Callback(hObject, eventdata, handles)
% hObject    handle to startgroup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

lead group

% --- Executes on button press in startmacaque.
function startmacaque_Callback(hObject, eventdata, handles)
% hObject    handle to startmacaque (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

lead macaque


% --- Executes on button press in startanatomy.
function startanatomy_Callback(hObject, eventdata, handles)
% hObject    handle to startanatomy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

lead anatomy
