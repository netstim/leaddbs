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

% Last Modified by GUIDE v2.5 29-Dec-2018 11:48:39

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

% list presets
earoot=ea_getearoot;
presf=[earoot,'ext_libs',filesep,'ANTs',filesep,'presets',filesep];
ndir=dir([presf,'ea_antspreset_*.m']);
namecell = cell(size(ndir));
funcell = cell(size(ndir));
for n=1:length(ndir)
    [~,funame,~]=fileparts(ndir(n).name);
    namecell{n}=eval([funame,'(''query'')']);
    funcell{n}=funame;
end
[namecell, sortix]=sort(namecell);
funcell = funcell(sortix);
setappdata(handles.pcpopup,'funcell',funcell);
set(handles.pcpopup,'String',namecell);

% select last selection
prefs=ea_prefs('');

% preset
[~,ix]=ismember(prefs.machine.normsettings.ants_preset,getappdata(handles.pcpopup,'funcell'));
if ix % if has prior selection
    set(handles.pcpopup,'Value',ix);
end

% metric
[~,ix]=ismember(prefs.machine.normsettings.ants_metric,get(handles.metric,'String'));
if ix % if has prior selection
    set(handles.metric,'Value',ix);
end

% strategy
[~,ix]=ismember(prefs.machine.normsettings.ants_strategy,get(handles.strategy,'String'));
if ix % if has prior selection
    set(handles.strategy,'Value',ix);
end

set(handles.includefa,'Value',prefs.machine.normsettings.ants_usefa);

if ischar(prefs.machine.normsettings.ants_numcores)
    set(handles.restrcores,'Value',1);
    set(handles.numcores,'String',prefs.machine.normsettings.ants_numcores)
else
    set(handles.restrcores,'Value',0);
    set(handles.numcores,'enable','off');
end

set(handles.scrf,'Value',prefs.machine.normsettings.ants_scrf);

set(handles.skullstripped,'Value',prefs.machine.normsettings.ants_skullstripped);

set(handles.usepreexisting,'Value',prefs.machine.normsettings.ants_usepreexisting);


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
normsettings.ants_preset=getappdata(handles.pcpopup,'funcell');
normsettings.ants_preset=normsettings.ants_preset{get(handles.pcpopup,'Value')};
normsettings.ants_scrf=get(handles.scrf,'Value');

normsettings.ants_metric=get(handles.metric,'String');
normsettings.ants_metric=normsettings.ants_metric{get(handles.metric,'Value')};

normsettings.ants_strategy=get(handles.strategy,'String');
normsettings.ants_strategy=normsettings.ants_strategy{get(handles.strategy,'Value')};

normsettings.ants_usefa=get(handles.includefa,'Value');

normsettings.ants_skullstripped=get(handles.skullstripped,'Value');

normsettings.ants_usepreexisting=get(handles.usepreexisting,'Value');


if get(handles.restrcores,'Value')
    normsettings.ants_numcores=get(handles.numcores,'String');
else
    normsettings.ants_numcores=0;
end

ea_setprefs('normsettings',normsettings);

delete(handles.setfig);


% --- Executes on button press in scrf.
function scrf_Callback(hObject, eventdata, handles)
% hObject    handle to scrf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of scrf

% --- Executes on selection change in metric.
function metric_Callback(hObject, eventdata, handles)
% hObject    handle to metric (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns metric contents as cell array
%        contents{get(hObject,'Value')} returns selected item from metric


% --- Executes during object creation, after setting all properties.
function metric_CreateFcn(hObject, eventdata, handles)
% hObject    handle to metric (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in strategy.
function strategy_Callback(hObject, eventdata, handles)
% hObject    handle to strategy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns strategy contents as cell array
%        contents{get(hObject,'Value')} returns selected item from strategy


% --- Executes during object creation, after setting all properties.
function strategy_CreateFcn(hObject, eventdata, handles)
% hObject    handle to strategy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function numcores_Callback(hObject, eventdata, handles)
% hObject    handle to numcores (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of numcores as text
%        str2double(get(hObject,'String')) returns contents of numcores as a double


% --- Executes during object creation, after setting all properties.
function numcores_CreateFcn(hObject, eventdata, handles)
% hObject    handle to numcores (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in restrcores.
function restrcores_Callback(hObject, eventdata, handles)
% hObject    handle to restrcores (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of restrcores
if get(handles.restrcores,'Value')
    set(handles.numcores,'enable','on');
else
    set(handles.numcores,'enable','off');
end

% --- Executes on button press in includefa.
function includefa_Callback(hObject, eventdata, handles)
% hObject    handle to includefa (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of includefa


% --- Executes on button press in skullstripped.
function skullstripped_Callback(hObject, eventdata, handles)
% hObject    handle to skullstripped (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of skullstripped

% --- Executes on selection change in usepreexisting.
function usepreexisting_Callback(hObject, eventdata, handles)
% hObject    handle to usepreexisting (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns usepreexisting contents as cell array
%        contents{get(hObject,'Value')} returns selected item from usepreexisting


% --- Executes during object creation, after setting all properties.
function usepreexisting_CreateFcn(hObject, eventdata, handles)
% hObject    handle to usepreexisting (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
