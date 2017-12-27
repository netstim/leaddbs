function varargout = ea_presults_horn2017(varargin)
% EA_PRESULTS_HORN2017 MATLAB code for ea_presults_horn2017.fig
%      EA_PRESULTS_HORN2017, by itself, creates a new EA_PRESULTS_HORN2017 or raises the existing
%      singleton*.
%
%      H = EA_PRESULTS_HORN2017 returns the handle to a new EA_PRESULTS_HORN2017 or the handle to
%      the existing singleton*.
%
%      EA_PRESULTS_HORN2017('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in EA_PRESULTS_HORN2017.M with the given input arguments.
%
%      EA_PRESULTS_HORN2017('Property','Value',...) creates a new EA_PRESULTS_HORN2017 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ea_presults_horn2017_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ea_presults_horn2017_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ea_presults_horn2017

% Last Modified by GUIDE v2.5 06-Dec-2017 17:37:53

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ea_presults_horn2017_OpeningFcn, ...
                   'gui_OutputFcn',  @ea_presults_horn2017_OutputFcn, ...
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


% --- Executes just before ea_presults_horn2017 is made visible.
function ea_presults_horn2017_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ea_presults_horn2017 (see VARARGIN)

% Choose default command line output for ea_presults_horn2017
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

cfg=varargin{1};
set(0,'CurrentFigure',handles.predfigure);
set(handles.predfigure,'CurrentAxes',handles.fMRImodel);
if isfield(cfg,'fMRI')
    imshow(cfg.fMRI.model);
end
axis off
set(handles.predfigure,'CurrentAxes',handles.fMRIpat);
if isfield(cfg,'fMRI')
    imshow(cfg.fMRI.vta);
end
axis off
set(handles.predfigure,'CurrentAxes',handles.dMRImodel);
if isfield(cfg,'dMRI')
    imshow(cfg.dMRI.model);
end
axis off
set(handles.predfigure,'CurrentAxes',handles.dMRIpat);
if isfield(cfg,'dMRI')
    imshow(cfg.dMRI.vta);
end
axis off


set(handles.stimlabel,'String',cfg.stim.label);
set(handles.patname,'String',cfg.stim.patientname);
set(handles.vtamodel,'String',cfg.stim.model);

set(handles.actCR,'String',num2str(find(cfg.stim.activecontacts{1})));
set(handles.actCL,'String',num2str(find(cfg.stim.activecontacts{2})));
set(handles.ampsR,'String',num2str((cfg.stim.amplitude{1}(logical(cfg.stim.amplitude{1})))));
set(handles.ampsL,'String',num2str((cfg.stim.amplitude{2}(logical(cfg.stim.amplitude{2})))));

set(handles.averror,'String',['+/- ',sprintf('%.02f',cfg.res.updrs3err),' %']);
set(handles.updrs3hat,'String',[sprintf('%.02f',cfg.res.updrs3imp),' %']);

set(handles.predfigure,'Name',['Predictions for ',cfg.stim.patientname,' (stimulation: ',cfg.stim.label,')']);

% UIWAIT makes ea_presults_horn2017 wait for user response (see UIRESUME)

% uiwait(handles.predfigure);


% --- Outputs from this function are returned to the command line.
function varargout = ea_presults_horn2017_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
