function varargout = ea_normsettings_maget_segment(varargin)
% EA_NORMSETTINGS_MAGET_SEGMENT MATLAB code for ea_normsettings_maget_segment.fig
%      EA_NORMSETTINGS_MAGET_SEGMENT, by itself, creates a new EA_NORMSETTINGS_MAGET_SEGMENT or raises the existing
%      singleton*.
%
%      H = EA_NORMSETTINGS_MAGET_SEGMENT returns the handle to a new EA_NORMSETTINGS_MAGET_SEGMENT or the handle to
%      the existing singleton*.
%
%      EA_NORMSETTINGS_MAGET_SEGMENT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in EA_NORMSETTINGS_MAGET_SEGMENT.M with the given input arguments.
%
%      EA_NORMSETTINGS_MAGET_SEGMENT('Property','Value',...) creates a new EA_NORMSETTINGS_MAGET_SEGMENT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ea_normsettings_maget_segment_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ea_normsettings_maget_segment_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ea_normsettings_maget_segment

% Last Modified by GUIDE v2.5 17-Aug-2016 14:47:46

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ea_normsettings_maget_segment_OpeningFcn, ...
                   'gui_OutputFcn',  @ea_normsettings_maget_segment_OutputFcn, ...
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


% --- Executes just before ea_normsettings_maget_segment is made visible.
function ea_normsettings_maget_segment_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ea_normsettings_maget_segment (see VARARGIN)

% Choose default command line output for ea_normsettings_maget_segment
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

set(handles.setfig,'Name','Normalization Settings');
set(handles.titletext,'String','MAGeT Normalization');


prefs=ea_prefs('');
[~,ix]=ismember(prefs.machine.normsettings.maget_peerset,get(handles.pcpopup,'String'));
set(handles.pcpopup,'Value',ix);

atlcells=get(handles.atlaspopup,'String');
if ~iscell(atlcells)
    atlcells={atlcells};
end
[~,ix]=ismember(prefs.machine.normsettings.maget_atlasset,atlcells);
set(handles.atlaspopup,'Value',ix);




% UIWAIT makes ea_normsettings_maget_segment wait for user response (see UIRESUME)


% --- Outputs from this function are returned to the command line.
function varargout = ea_normsettings_maget_segment_OutputFcn(hObject, eventdata, handles) 
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
normsettings.maget_peerset=get(handles.pcpopup,'String');
normsettings.maget_peerset=normsettings.maget_peerset{get(handles.pcpopup,'Value')};
normsettings.maget_peersetcell=getappdata(handles.pcpopup,'peersetcell');
normsettings.maget_atlasset=get(handles.atlaspopup,'String');

if iscell(normsettings.maget_atlasset)
normsettings.maget_atlasset=normsettings.maget_atlasset{get(handles.atlaspopup,'Value')};
end

ea_setprefs('normsettings',normsettings);
delete(handles.setfig);


% --- Executes on selection change in atlaspopup.
function atlaspopup_Callback(hObject, eventdata, handles)
% hObject    handle to atlaspopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns atlaspopup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from atlaspopup


% --- Executes during object creation, after setting all properties.
function atlaspopup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to atlaspopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
