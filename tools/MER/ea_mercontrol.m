function varargout = ea_mercontrol(varargin)
%EA_MERCONTROL MATLAB code file for ea_mercontrol.fig
%      EA_MERCONTROL, by itself, creates a new EA_MERCONTROL or raises the existing
%      singleton*.
%
%      H = EA_MERCONTROL returns the handle to a new EA_MERCONTROL or the handle to
%      the existing singleton*.
%
%      EA_MERCONTROL('Property','Value',...) creates a new EA_MERCONTROL using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to ea_mercontrol_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      EA_MERCONTROL('CALLBACK') and EA_MERCONTROL('CALLBACK',hObject,...) call the
%      local function named CALLBACK in EA_MERCONTROL.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ea_mercontrol

% Last Modified by GUIDE v2.5 07-Jul-2017 14:57:26

% __________________________________________________________________________________
% Copyright (C) 2017 University of Pittsburgh, Brain Modulation Lab
%
% Ari Kappel

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ea_mercontrol_OpeningFcn, ...
                   'gui_OutputFcn',  @ea_mercontrol_OutputFcn, ...
                   'gui_LayoutFcn',  [], ...
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





%%%%%%%%%%%%%%%%%% GUIDE-generated functions  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes just before ea_mercontrol is made visible.
function ea_mercontrol_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)
clc
set(hObject, 'Name', 'MER Control',...
    'CloseRequestFcn', @closefunction, 'KeyPressFcn', @ea_keypress)
% Choose default command line output for ea_mercontrol
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ea_mercontrol wait for user response (see UIRESUME)
% uiwait(handles.mercontrolfig);

resultfig=varargin{1};
options=varargin{2};
setappdata(hObject, 'resultfig', resultfig);
setappdata(hObject, 'options', options);  % Keep a copy of options
setappdata(resultfig, 'merstruct', options.prefs.mer);  %.length [24 mm]; .offset [2 mm]; .tract_info

ea_gui_generate(handles);  % Generate UI elements for each tract.
ea_resultfig_clear(handles);  % Delete markers and mer trajectories.
ea_data_clear(handles);  % Clears trajectories, depths, etc in merstruct.
ea_data_setdefault(handles);  % Sets merstruct and mermarkers to default values.
ea_resultfig_update(handles);  % Plot markers and trajs from data.
ea_mercontrol_updateall(handles);  % Set GUI from data.


% --- Outputs from this function are returned to the command line.
function varargout = ea_mercontrol_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in pushtrackmer_left.
function pushtrackmer_left_Callback(hObject, eventdata, handles)
% hObject    handle to pushtrackmer_left (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns pushtrackmer_left contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pushtrackmer_left
sidestr='left';
resultfig=getappdata(handles.mercontrolfig,'resultfig');
options=getappdata(handles.mercontrolfig,'options');
ea_trackmer(sidestr,resultfig,options)


% --- Executes during object creation, after setting all properties.
function pushtrackmer_left_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pushtrackmer_left (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in pushimportrec_left.
function pushimportrec_left_Callback(hObject, eventdata, handles)
% hObject    handle to pushimportrec_left (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns pushimportrec_left contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pushimportrec_left
sidestr='left'; side=2;
resultfig=getappdata(handles.mercontrolfig,'resultfig');
options=getappdata(handles.mercontrolfig,'options');
[Depth,Thresh,Tracts] = ea_readmer(fullfile(options.uipatdirs{1},['Rec_',sidestr,'.mat']));
merlfp.(sidestr).Depth = Depth;
merlfp.(sidestr).Tracts = Tracts;
merlfp.(sidestr).Thresh = Thresh;

save(fullfile(options.uipatdirs{1},'ea_recordings.mat'),'merlfp')
disp('Process DONE')

% --- Executes during object creation, after setting all properties.
function pushimportrec_left_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pushimportrec_left (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmedial_left.
function popupmedial_left_Callback(hObject, eventdata, handles)
% hObject    handle to popupmedial_left (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmedial_left contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmedial_left


% --- Executes during object creation, after setting all properties.
function popupmedial_left_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmedial_left (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in pushtrackmer_right.
function pushtrackmer_right_Callback(hObject, eventdata, handles)
% hObject    handle to pushtrackmer_right (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns pushtrackmer_right contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pushtrackmer_right


% --- Executes during object creation, after setting all properties.
function pushtrackmer_right_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pushtrackmer_right (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in pushimportrec_right.
function pushimportrec_right_Callback(hObject, eventdata, handles)
% hObject    handle to pushimportrec_right (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns pushimportrec_right contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pushimportrec_right


% --- Executes during object creation, after setting all properties.
function pushimportrec_right_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pushimportrec_right (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmedial_right.
function popupmedial_right_Callback(hObject, eventdata, handles)
% hObject    handle to popupmedial_right (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmedial_right contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmedial_right


% --- Executes during object creation, after setting all properties.
function popupmedial_right_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmedial_right (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on selection change in popupimplantedtract_left.
function popupimplantedtract_left_Callback(hObject, eventdata, handles)
popupimplantedtract_helper(hObject, eventdata, handles);


% --- Executes during object creation, after setting all properties.
function popupimplantedtract_left_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupimplantedtract_left (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupimplantedtract_right.
function popupimplantedtract_right_Callback(hObject, eventdata, handles)
popupimplantedtract_helper(hObject, eventdata, handles);


% --- Executes during object creation, after setting all properties.
function popupimplantedtract_right_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupimplantedtract_right (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function editimplanteddepth_left_Callback(hObject, eventdata, handles)
editimplanteddepth_helper(hObject, eventdata, handles);


% --- Executes during object creation, after setting all properties.
function editimplanteddepth_left_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editimplanteddepth_left (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function editimplanteddepth_right_Callback(hObject, eventdata, handles)
editimplanteddepth_helper(hObject, eventdata, handles);


% --- Executes during object creation, after setting all properties.
function editimplanteddepth_right_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editimplanteddepth_right (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmermarkers_right.
function popupmermarkers_right_Callback(hObject, eventdata, handles)
% hObject    handle to popupmermarkers_right (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmermarkers_right contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmermarkers_right


% --- Executes during object creation, after setting all properties.
function popupmermarkers_right_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmermarkers_right (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmermarkers_left.
function popupmermarkers_left_Callback(hObject, eventdata, handles)
% hObject    handle to popupmermarkers_left (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmermarkers_left contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmermarkers_left


% --- Executes during object creation, after setting all properties.
function popupmermarkers_left_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmermarkers_left (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in clearall.
function clearall_Callback(hObject, eventdata, handles)
% hObject    handle to clearall (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ea_resultfig_clear(handles);  % Delete markers and mer trajectories.
ea_data_clear(handles);  % Clears trajectories, depths, etc in merstruct.
ea_mercontrol_updateall(handles);  % Set GUI from data.


% --- Executes during object creation, after setting all properties.
function clearall_CreateFcn(hObject, eventdata, handles)
% hObject    handle to clearall (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
root = fileparts(which('ea_coregmr.m'));
background = imread([root,'/icons/delete.png']);
set(hObject,'CData',background);


% --- Executes on button press in setdefault.
function setdefault_Callback(hObject, eventdata, handles)
% hObject    handle to setdefault (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ea_resultfig_clear(handles);  % Delete markers and mer trajectories.
ea_data_clear(handles);  % Clears trajectories, depths, etc in merstruct.
ea_data_setdefault(handles);  % Sets merstruct and mermarkers to default values.
ea_resultfig_update(handles);  % Plot markers and trajs from data.
ea_mercontrol_updateall(handles);  % Set GUI from data.


% --- Executes during object creation, after setting all properties.
function setdefault_CreateFcn(hObject, eventdata, handles)
% hObject    handle to setdefault (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
root = fileparts(which('ea_coregmr.m'));
background = imread([root,'/icons/checkmark.png']);
set(hObject,'CData',background)


% --- Executes on button press in undomarker.
function undomarker_Callback(hObject, eventdata, handles)
% hObject    handle to undomarker (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
resultfig = getappdata(handles.mercontrolfig, 'resultfig');
mermarkers = getappdata(resultfig, 'mermarkers');
markerhistory = getappdata(resultfig, 'markerhistory');
markerhistory(end+1) = mermarkers(end);
delete(markerhistory(end).handle);
markerhistory(end).handle = [];
delete(markerhistory(end).tag.handle);
markerhistory(end).tag.handle = [];
mermarkers(end) = [];
setappdata(resultfig, 'mermarkers', mermarkers);
setappdata(resultfig, 'markerhistory', markerhistory);
ea_resultfig_updatemarkers(handles);
ea_mercontrol_updatemarkers(handles);


% --- Executes during object creation, after setting all properties.
function undomarker_CreateFcn(hObject, eventdata, handles)
% hObject    handle to undomarker (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
root = fileparts(which('ea_coregmr.m'));
background = imread([root,'/icons/undo.png']);
set(hObject,'CData',background,'Value',0)


% --- Executes on button press in redomarker.
function redomarker_Callback(hObject, eventdata, handles)
% hObject    handle to redomarker (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
resultfig = getappdata(handles.mercontrolfig, 'resultfig');
mermarkers = getappdata(resultfig, 'mermarkers');
markerhistory = getappdata(resultfig, 'markerhistory');
mermarkers(end + 1) = markerhistory(end);
markerhistory(end) = [];
setappdata(resultfig, 'mermarkers', mermarkers);
setappdata(resultfig, 'markerhistory', markerhistory);
ea_resultfig_updatemarkers(handles);
ea_mercontrol_updatemarkers(handles);


% --- Executes during object creation, after setting all properties.
function redomarker_CreateFcn(hObject, eventdata, handles)
% hObject    handle to redomarker (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
root = fileparts(which('ea_coregmr.m'));
background = imread([root,'/icons/redo.png']);
set(hObject,'CData',background,'Value',0)


% --- Executes on button press in togglemarkertags.
function togglemarkertags_Callback(hObject, eventdata, handles)
% hObject    handle to togglemarkertags (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglemarkertags
resultfig = getappdata(handles.mercontrolfig,'resultfig');
merstruct = getappdata(resultfig, 'merstruct');
if get(hObject,'Value')
    merstruct.tag.visible = 'on';
else
    merstruct.tag.visible = 'off';
end
setappdata(resultfig, 'merstruct', merstruct);
ea_resultfig_updatemarkers(handles);


% --- Executes during object creation, after setting all properties.
function togglemarkertags_CreateFcn(hObject, eventdata, handles)
% hObject    handle to togglemarkertags (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
root = fileparts(which('ea_coregmr.m'));
background = imread([root,'/icons/text.png']);
set(hObject,'CData',background,'Value',0)


% --- Executes on button press in importmarkers.
function importmarkers_Callback(hObject, eventdata, handles)
% hObject    handle to importmarkers (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
resultfig = getappdata(handles.mercontrolfig, 'resultfig');
options = getappdata(handles.mercontrolfig, 'options');

try 
    load([options.uipatdirs{1},filesep,'ea_mermarkers.mat'],...
        'mermarkers', 'mertoggles');
catch
    [files] = fdir(options.uipatdirs{1},'.mat');
    if ~isempty(cell2mat(strfind({files(:).name},'ea_mermarkers')))
        %[s,b] = listdlg('PromptString','Select a file:','SelectionMode','single',...
        %   'ListString',{files(~cellfun(@isempty,strfind({files(:).name},'ea_mermarkers'))).name})
    end
    ea_error(['could not find ea_mermarkers.mat for ', options.patientname])
end
ea_resultfig_clear_markers;
setappdata(resultfig, 'mermarkers', mermarkers);
setappdata(handles.mercontrolfig, 'mertoggles', mertoggles);
ea_resultfig_updatemarkers(handles);
ea_mercontrol_updatemarkers(handles);
ea_mercontrol_updatetoggles(handles);


% --- Executes during object creation, after setting all properties.
function importmarkers_CreateFcn(hObject, eventdata, handles)
% hObject    handle to importmarkers (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
root = fileparts(which('ea_coregmr.m'));
background = imread([root,'/icons/import.png']);
set(hObject,'CData',background,'Value',0)


% --- Executes on button press in exportmarkers.
function exportmarkers_Callback(hObject, eventdata, handles)
% hObject    handle to exportmarkers (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
resultfig = getappdata(handles.mercontrolfig, 'resultfig');
options = getappdata(handles.mercontrolfig, 'options');
mermarkers = getappdata(resultfig, 'mermarkers');
mertoggles = getappdata(handles.mercontrolfig, 'mertoggles'); %#ok<NASGU>

for marker_ix = 1:length(mermarkers)
    % Empty the handle to make sure it gets redrawn on import.
    mermarkers(marker_ix).handle = [];
    if isfield(mermarkers(marker_ix).tag, 'handle')
        mermarkers(marker_ix).tag.handle = [];
    end
end

if exist([options.uipatdirs{1},filesep,'ea_mermarkers.mat'],'file')
    overwrite = ea_questdlg({['ea_mermarkers.mat found in ' options.patientname ' directory.'];...
        'Would you like to overwrite this file?'}, options.patientname);
end
if ~exist([options.uipatdirs{1},filesep,'ea_mermarkers.mat'],'file') || strcmpi(overwrite,'Yes')
    disp(['Saving: ', options.uipatdirs{1}, filesep, 'ea_mermarkers.mat']);
    save([options.uipatdirs{1}, filesep, 'ea_mermarkers.mat'],...
        'mermarkers','mertoggles');
    disp('DONE');
end


% --- Executes during object creation, after setting all properties.
function exportmarkers_CreateFcn(hObject, eventdata, handles)
% hObject    handle to exportmarkers (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
root = fileparts(which('ea_coregmr.m'));
background = imread([root,'/icons/export.png']);
set(hObject,'CData',background,'Value',0)


% --- Executes on button press in pushbutton_clear_checks.
function pushbutton_clear_checks_Callback(hObject, eventdata, handles)
mertoggles = getappdata(handles.mercontrolfig, 'mertoggles');
mertoggles.keycontrol = zeros(2, length(merstruct.tract_info));
setappdata(handles.mercontrolfig, 'mertoggles', mertoggles);
ea_mercontrol_updatetoggles(handles);
















%%%%%%%%%%%%%%%%% Custom functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ea_keypress(handles, event) %#ok<INUSL>
% this is the main keypress function for the mercontrolfig.
% Add event listeners here.
if ismember('shift', event.Modifier)
end
fprintf('Keypress ignored. Click on resultfig first.\n')


function closefunction(hObject, event) %#ok<INUSD>
handles = getappdata(hObject,'UsedByGUIData_m');
ea_resultfig_clear(handles);
ea_data_clear(handles);
delete(hObject)


function ea_gui_generate(handles)
% ea_gui_generate(handles)
% For each side, for each tract_inf, generates
% -display toggle button
% -position/depth edit box
% -keyboard control toggle

merstruct = getappdata(getappdata(handles.mercontrolfig, 'resultfig'), 'merstruct');
% Generate UI elements
width = handles.mainuipanel.InnerPosition(3);
popup_string = cat(2, 'Set implanted...', {merstruct.tract_info.label});
for side_str = {'left', 'right'}
    % In mainuipanel...
    
    % ...Add options to popupimplantedtract_left _right
    set(handles.(['popupimplantedtract_' side_str{:}]),...
        'String', popup_string);
    
    % ...Add togglebuttons next to texttoggle_left and _right (with callbacks)
    anchor_pos = get(handles.(['texttoggle_' side_str{:}]), 'Position');
    w_per = (width - anchor_pos(3) - anchor_pos(1)) / length(merstruct.tract_info);
    for tract_ix = 1:length(merstruct.tract_info)
        pos_str = merstruct.tract_info(tract_ix).label;
        new_pos = [anchor_pos(1) + anchor_pos(3) + (tract_ix-1)*w_per...
            anchor_pos(2)-5 w_per 26];
        uicontrol(handles.mainuipanel,...
            'Style', 'togglebutton',...
            'Tag', ['toggle' pos_str '_' side_str{:}],...
            'Value', 1,...
            'String', pos_str,...
            'Position', new_pos,...
            'Units', 'pixels',...
            'BackgroundColor', [1 1 1],...
            'Value', 1,...
            'Callback', @ea_updatetogglestate);
    end
    
    % ...Add edit boxes next to textpos_left and _right (with callbacks)
    anchor_pos = get(handles.(['textpos_' side_str{:}]), 'Position');
    w_per = (width - anchor_pos(3) - anchor_pos(1)) / length(merstruct.tract_info);
    for tract_ix = 1:length(merstruct.tract_info)
        pos_str = merstruct.tract_info(tract_ix).label;
        new_pos = [anchor_pos(1) + anchor_pos(3) + (tract_ix-1)*w_per...
            anchor_pos(2)-5 w_per 26];
        uicontrol(handles.mainuipanel,...
            'Style', 'edit',...
            'Tag', ['pos' pos_str '_' side_str{:}],...
            'String', '0',...
            'Position', new_pos,...
            'Units', 'pixels',...
            'BackgroundColor', [1 1 1],...
            'Callback', @ea_updatepostext);
    end
    
    % In keyboardcontrolgroup...
    % ...Add radiobuttons next to keytext_left and _right inside
    anchor_pos = get(handles.(['textkey_' side_str{:}]), 'Position');
    w_per = (width - anchor_pos(3) - anchor_pos(1)) / length(merstruct.tract_info);
    for tract_ix = 1:length(merstruct.tract_info)
        pos_str = merstruct.tract_info(tract_ix).label;
        new_pos = [anchor_pos(1) + anchor_pos(3) + (tract_ix-1)*w_per...
            anchor_pos(2) w_per anchor_pos(4)];
        uicontrol(handles.keyboardcontrolgroup,...
            'Style', 'checkbox',...
            'Tag', ['keycheck' pos_str '_' side_str{:}],...
            'String', pos_str,...
            'Position', new_pos,...
            'Units', 'pixels',...
            'BackgroundColor', [1 1 1],...
            'Callback', @ea_keycheck);
    end
end


function ea_keycheck(hObject, eventdata) %#ok<INUSD>
[~, sid, track] = ea_detsidestr(hObject.Tag);
% Update the value
handles = guidata(hObject);
resultfig = getappdata(handles.mercontrolfig, 'resultfig');
merstruct = getappdata(resultfig, 'merstruct');
bTrack = strcmpi({merstruct.tract_info.label}, track);
mertoggles = getappdata(handles.mercontrolfig, 'mertoggles');
mertoggles.keycontrol(sid, bTrack) = get(hObject,'Value') == get(hObject,'Max');
setappdata(handles.mercontrolfig, 'mertoggles', mertoggles);


function ea_data_clear(handles, varargin)
% ea_clear_data(handles)
% Resets trajectories, depths, etc.
ea_data_clear_traj(handles, varargin{:});
ea_data_clear_markers(handles);


function ea_data_clear_traj(handles, sidestr)
resultfig = getappdata(handles.mercontrolfig, 'resultfig');
merstruct = getappdata(resultfig, 'merstruct');
% if ~exist('sidestr','var')
%     sidestr = 'both';
% end
% [side_strs, clear_sides, ~] = ea_detsidestr(sidestr);
if isfield(merstruct, 'currentmer')
    merstruct = rmfield(merstruct, 'currentmer');
end
setappdata(resultfig, 'merstruct', merstruct);


function ea_data_clear_markers(handles)
resultfig = getappdata(handles.mercontrolfig, 'resultfig');
mermarkers = struct('side', {}, 'tract', {}, 'depth', {},...
    'markertype', {}, 'session', {}, 'dat', {}, 'tag', {},...
    'handle', {}, 'notes', {}, 'coords_mm', {});
setappdata(resultfig, 'mermarkers', mermarkers);
setappdata(resultfig, 'markerhistory', mermarkers);


function ea_resultfig_clear(handles, sidestr)
% ea_resultfig_clear(handles)
% Delete markers and plot empty 3dplot trajs.
if ~exist('sidestr','var')
    sidestr = 'both';
end
ea_resultfig_clear_traj(handles, sidestr);
ea_resultfig_clear_markers(handles, sidestr);


function ea_resultfig_clear_traj(handles, sidestr)
resultfig = getappdata(handles.mercontrolfig, 'resultfig');
if ~exist('sidestr','var')
    sidestr = 'both';
end
[side_strs, clear_sides, ~] = ea_detsidestr(sidestr);

% Get a list of tags for all the plotted items in resultfig.
curr_fig_items = resultfig.CurrentAxes.Children;
tag_array = cell(1, length(curr_fig_items));
for cfi = 1:length(curr_fig_items)
    tag_array{cfi} = curr_fig_items(cfi).Tag;
end

merhandles = getappdata(resultfig,'merhandles');
merstruct = getappdata(resultfig,'merstruct');
set(0, 'CurrentFigure', resultfig); hold on;
for sid = clear_sides
    sidestr = side_strs{sid};
    for pos_ix = 1:length(merstruct.tract_info)
        pos_str = merstruct.tract_info(pos_ix).label;
        try
            if isstruct(merhandles) && isfield(merhandles, pos_str) &&...
                    length(merhandles.(pos_str)) >= sid &&...
                    ~isempty(merhandles.(pos_str){sid})
                delete(merhandles.(pos_str){sid})
            end
            resid_items = curr_fig_items(strcmpi(tag_array, [pos_str '_' sidestr]));
            if ~isempty(resid_items)
                delete(resid_items);
            end
        catch
            fprintf('TODO: Warning message about deleting handles.\n')
        end
        % Populate merhandles. Because we aren't plotting anything, the
        % object handle will remain empty. color and tag are ignored.
        merhandles.(pos_str){sid} = plot3([], [], [],...
            'color', merstruct.tract_info(pos_ix).color, 'linew', 5,...
            'tag', [pos_str '_' sidestr]);
    end
end
setappdata(resultfig, 'merhandles', merhandles);


function ea_resultfig_clear_markers(handles, sidestr)
% ea_resultfig_clear_markers(handles, sidestr)
if ~exist('sidestr','var')
    sidestr = 'both';
end
[side_strs, ~, ~] = ea_detsidestr(sidestr);
resultfig = getappdata(handles.mercontrolfig, 'resultfig');
mermarkers = getappdata(resultfig,'mermarkers');
for marker_ix = 1:length(mermarkers)
    if any(strcmpi(mermarkers(marker_ix).side, side_strs))
        try delete(mermarkers(marker_ix).handle); end
        try delete(mermarkers(marker_ix).tag.handle); end
        mermarkers(marker_ix).handle = [];
        if isfield(mermarkers(marker_ix).tag, 'handle')
            mermarkers(marker_ix).tag.handle = [];
        end
    end
end
setappdata(resultfig,'mermarkers',mermarkers)


function ea_data_setdefault(handles)
% ea_data_setdefault(handles)
% Sets data model to default values.
ea_data_setdefault_trajs(handles);
ea_data_setdefault_markers(handles);
ea_data_setdefault_mertoggles(handles);


function ea_data_setdefault_trajs(handles)
% ea_data_setdefault_trajs(handles)
resultfig = getappdata(handles.mercontrolfig, 'resultfig');
options = getappdata(handles.mercontrolfig, 'options');
merstruct = getappdata(resultfig, 'merstruct');
[~, ~, dbs_contacts] = ea_load_reconstruction(options);
merstruct.dbs_contacts_mm = ea_resolvecoords(dbs_contacts, options);
merstruct.implant_depth = [0 0];
merstruct.implant_idx = merstruct.defaulttract * [1 1];
[~, side_ids, ~] = ea_detsidestr('both');
for sid = side_ids
    for pos_ix = 1:length(merstruct.tract_info)
        pos_str = merstruct.tract_info(pos_ix).label;
        merstruct.currentmer.(pos_str).dist(sid) = 0;
    end
end
setappdata(resultfig, 'merstruct', merstruct);
ea_mercontrol_updatetrajectories(handles);  % Update the trajectories stored in merstruct.


function ea_data_setdefault_markers(handles)
% ea_data_setdefault_markers(handles)
% Simply clears markers. See ea_data_clear_markers.
ea_data_clear_markers(handles);


function ea_data_setdefault_mertoggles(handles)
resultfig = getappdata(handles.mercontrolfig, 'resultfig');
merstruct = getappdata(resultfig, 'merstruct');
mertoggles.keycontrol = zeros(2, length(merstruct.tract_info));
mertoggles.togglestates = ones(2, length(merstruct.tract_info));
setappdata(handles.mercontrolfig, 'mertoggles', mertoggles);


function ea_resultfig_update(handles)
% ea_resultfig_update(handles)
% Updates trajectory plot and markers plot from data.
ea_resultfig_updatetrajectories(handles);
ea_resultfig_updatemarkers(handles);


function ea_mercontrol_updateall(handles)
% ea_mercontrol_updateall(handles)
% Update GUI elements from data structs.
ea_mercontrol_updateimplanted(handles);
ea_mercontrol_updatetoggles(handles);
ea_mercontrol_updatemarkers(handles);


function ea_mercontrol_updatetoggles(handles, side_str)
if ~exist('side_str', 'var')
    side_str = 'both';
end
[side_strs, side_ids, ~] = ea_detsidestr(side_str);
resultfig = getappdata(handles.mercontrolfig, 'resultfig');
mertoggles = getappdata(handles.mercontrolfig, 'mertoggles');
merstruct = getappdata(resultfig, 'merstruct');
ui_tags = {handles.mainuipanel.Children.Tag};
kcg_tags = {handles.keyboardcontrolgroup.Children.Tag};
for sid = side_ids
    side_str = side_strs{sid};
    for tract_ix = 1:length(merstruct.tract_info)
        tract_info = merstruct.tract_info(tract_ix);
        val = mertoggles.togglestates(sid, tract_ix);
        set(handles.mainuipanel.Children(strcmpi(ui_tags,...
            ['toggle' tract_info.label '_' side_str])),...
            'Value', val);
        kh = handles.keyboardcontrolgroup.Children(...
            strcmpi(kcg_tags, ['keycheck' tract_info.label '_' side_str]));
        set(kh, 'Value', mertoggles.keycontrol(sid, tract_ix));
    end
end


function editimplanteddepth_helper(hObject, eventdata, handles)
dist = str2double(get(hObject, 'String'));
str_parts = strsplit(eventdata.Source.Tag, '_');
[sidestr, side_ix, ~] = ea_detsidestr(str_parts{2});
resultfig = getappdata(handles.mercontrolfig, 'resultfig');
merstruct = getappdata(resultfig, 'merstruct');
merstruct.implant_depth(side_ix) = dist;
for tract_ix = 1:length(merstruct.tract_info)
    pos_str = merstruct.tract_info(tract_ix).label;
    merstruct.currentmer.(pos_str).dist(side_ix) = dist;
end
setappdata(resultfig, 'merstruct', merstruct);
ea_mercontrol_updateimplanted(handles, sidestr{side_ix});
ea_resultfig_updatetrajectories(handles, sidestr{side_ix});


function popupimplantedtract_helper(hObject, eventdata, handles) %#ok<INUSL>
tract_ix = get(hObject, 'Value') - 1;
str_parts = strsplit(hObject.Tag, '_');
[~, sid, ~] = ea_detsidestr(str_parts{2});
resultfig = getappdata(handles.mercontrolfig, 'resultfig');
merstruct = getappdata(resultfig, 'merstruct');
merstruct.implant_idx(sid) = tract_ix;
setappdata(resultfig, 'merstruct', merstruct);
%TODO: I don't think the knowing the implanted trajectory actually matters.
% If it does, then...
% ea_mercontrol_updatetrajectories(handles, side_strs{sid});
% ea_resultfig_updatetrajectories();


function ea_updatepostext(hObject, eventData) %#ok<INUSD>
% Called when one of the position text boxes is edited.
handles = guidata(hObject);
resultfig = getappdata(handles.mercontrolfig, 'resultfig');
merstruct = getappdata(resultfig, 'merstruct');
[side_strs, sid, track] = ea_detsidestr(hObject.Tag);
merstruct.currentmer.(track).dist(sid) = str2double(get(hObject, 'String'));
setappdata(resultfig, 'merstruct', merstruct);
ea_mercontrol_updatetrajectories(handles, side_strs{sid});
ea_resultfig_updatetrajectories(handles);


function mertoggles = ea_updatetogglestate(hObject, eventdata) %#ok<INUSD>
% Called when one of the tract toggle buttons is pressed.
[~, sid, track] = ea_detsidestr(hObject.Tag);
handles = guidata(hObject);
merstruct = getappdata(getappdata(handles.mercontrolfig, 'resultfig'), 'merstruct');
bTrack = strcmpi({merstruct.tract_info.label}, track);
% Update the value in mertoggles
mertoggles = getappdata(handles.mercontrolfig, 'mertoggles');
mertoggles.togglestates(sid, bTrack) = hObject.Value;
setappdata(handles.mercontrolfig, 'mertoggles', mertoggles);
% Update the visual presentation
ea_resultfig_updatetrajectories(handles);


function [files,subs] = fdir(Path,ext)
% Returns contents of path as struct
% Option to keep only files with input extension
%
% Inputs:
%
% 
% Ari Kappel, 2017
%

if isempty(Path) || ~exist('Path','var')
    Path = pwd;
end
    
contents = dir(Path);
contents = contents(cellfun(@(x) isempty(regexp(x, '^\.', 'once')), {contents.name}));

if ~isempty({contents(~[contents.isdir]).name})
    files = contents(~[contents.isdir]);
end
if ~isempty({contents([contents.isdir]).name})
    subs = contents([contents.isdir]);
end

if nargin==2
    files = files(~cellfun(@isempty, strfind({files.name},ext)));
end
