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

% set lead path first time. to overwrite path run: lead path
if ~isdeployed && ~contains(path, fullfile(fileparts(mfilename),'helpers'))
    ea_setpath;
    if nargin>3
        switch varargin{1}
            case 'path'
                delete(hObject)
                return
        end
    end
end

ea_compat_data;
earoot=ea_getearoot;

ea_checkbuildspace;

ea_checkspm;

% check existence of directories
ea_checkleaddirs;

% check for commands first
if nargin == 4
    switch lower(varargin{1})
        case {'dbs', '-d', 'd'}
            lead_dbs;
            delete(handles.leadfigure)
            return
        case {'demo'}
            lead_dbs;
            lead_demo;
            delete(handles.leadfigure)
            return
        case {'group', '-g', 'g'}
            lead_group;
            delete(handles.leadfigure)
            return
        case {'groupconnectome', '-gc', 'gc'}
            lead_group_connectome;
            delete(handles.leadfigure)
            return
        case {'connectome', 'conn', '-c', 'c'}
            lead_connectome;
            delete(handles.leadfigure)
            return
        case {'mapper','connectomemapper', '-m', 'm'}
            lead_mapper;
            delete(handles.leadfigure)
            return
        case {'napper'}
            try % easter egg
                load('ea_napper.mat');
                sound(A,fs);
            end
        case {'or', '-o', 'o'}
            lead_or;
            delete(handles.leadfigure)
            return
        case {'anatomy', '-a', 'a'}
            lead_anatomy;
            delete(handles.leadfigure)
            return
        case {'predict', '-p', 'p'}
            lead_predict;
            delete(handles.leadfigure)
            return
        case {'import', '-i', 'i'}
            lead_import;
            delete(handles.leadfigure)
            return
        case {'version', 'ver', '-v', 'v'}
            disp(ea_getvsn('local'));
            delete(handles.leadfigure)
            return
        case 'dir'
            cd(ea_getearoot);
            delete(handles.leadfigure)
            return
        case 'path'
            ea_setpath;
            delete(handles.leadfigure)
            return
        case 'speak'
            fprintf('\n \n \n \n %s \n \n','L337-D8Z: "H3LLo 7o joO MY Phr13nd. l1V3 Lon9 4nd pRO5P3r."'); % yes, this indeed is an easter-egg.
            delete(handles.leadfigure)
            return
    end
elseif nargin == 5 && strcmp(varargin{1}, 'execute') % execute options specified in .json file
    set(hObject,'Visible','off'); drawnow;
    fid = fopen(varargin{2},'r');
    options = jsondecode(fread(fid,'*char')'); fclose(fid);
    ea_run('run',options);
    set(hObject,'Visible','on'); drawnow;
    close(hObject)
    delete(handles.leadfigure)
    return
elseif nargin > 5
    ea_command_line_run(varargin{:})
    delete(handles.leadfigure)
    return            
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

% Disable buttons for standalone app
% if isdeployed
%     set(handles.startconnectome,'Enable','off')
%     set(handles.startgroup,'Enable','off')
%     set(handles.startanatomy,'Enable','off')
% end

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
if isvalid(hObject) && nargout
    varargout{1} = hObject;
end

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

% --- Executes on button press in startanatomy.
function startanatomy_Callback(hObject, eventdata, handles)
% hObject    handle to startanatomy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

lead anatomy
