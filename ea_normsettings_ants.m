function varargout = ea_normsettings_ants(varargin)
% EA_NORMSETTINGS_ANTS MATLAB code for ea_normsettings_ants.fig
%      EA_NORMSETTINGS_ANTS, by itself, creates a new EA_NORMSETTINGS_ANTS or raises the existing
%      singleton*.
%
%      H = EA_NORMSETTINGS_ANTS returns the handle to a new EA_NORMSETTINGS_ANTS or the handle to
%      the existing singleton*.
%
%      EA_NORMSETTINGS_ANTS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in EA_NORMSETTINGS_ANTS.M with the given input arguments.
%
%      EA_NORMSETTINGS_ANTS('Property','Value',...) creates a new EA_NORMSETTINGS_ANTS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ea_normsettings_ants_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ea_normsettings_ants_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ea_normsettings_ants

% Last Modified by GUIDE v2.5 14-Apr-2017 19:07:39

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ea_normsettings_ants_OpeningFcn, ...
                   'gui_OutputFcn',  @ea_normsettings_ants_OutputFcn, ...
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


% --- Executes just before ea_normsettings_ants is made visible.
function ea_normsettings_ants_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ea_normsettings_ants (see VARARGIN)

% Choose default command line output for ea_normsettings_ants
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

set(handles.setfig,'Name','Normalization Settings');
set(handles.titletext,'String','ANTs Defaults');

prefs=ea_prefs('');
[~,ix]=ismember(prefs.machine.normsettings.ants_synmode,get(handles.pcpopup,'String'));
set(handles.pcpopup,'Value',ix);



% UIWAIT makes ea_normsettings_ants wait for user response (see UIRESUME)


% --- Outputs from this function are returned to the command line.
function varargout = ea_normsettings_ants_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure

% --- Executes on selection change in pcpopup.
function pcpopup_Callback(hObject, eventdata, handles)
% hObject    handle to pcpopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns pcpopup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pcpopup
selpc=get(hObject,'String');
selpc=selpc{get(hObject,'Value')};
if strcmp(selpc,'User-Defined')
   uipeerdir=ea_getpatients;
   setappdata(hObject,'peersetcell',uipeerdir);
end


% --- Executes during object creation, after setting all properties.
function pcpopup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pcpopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in savebutn.
function savebutn_Callback(hObject, eventdata, handles)
% hObject    handle to savebutn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
prefs=ea_prefs('');
normsettings=prefs.machine.normsettings;
normsettings.ants_synmode=get(handles.pcpopup,'String');
normsettings.ants_synmode=normsettings.ants_synmode{get(handles.pcpopup,'Value')};

ea_setprefs('normsettings',normsettings);

delete(handles.setfig);
