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

% Last Modified by GUIDE v2.5 28-Jul-2017 15:26:39

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

if nargin==4 % trajectory object supplied
    trajectory=varargin{1};
    resultfig=trajectory.plotFigureH;
    options=trajectory.options;
    setappdata(hObject,'trajectory',trajectory);
else % classical way of opening MER control figure
    resultfig=varargin{1};
    options=varargin{2};
end
setappdata(hObject, 'resultfig', resultfig);
setappdata(hObject, 'options', options);  % Keep a copy of options

% Get the MER State, set it to defaults if not supplied, and store it in appdata.
if exist('trajectory','var') && ~isempty(trajectory.merstruct)
    merstruct=trajectory.merstruct;
else
    merstruct = MERState();
    merstruct.setOptions(options);
    merstruct.clearData();
    if exist('trajectory','var')
        merstruct.Trajectory=trajectory;
    end
    merstruct.setDataToDefaults();
end

setappdata(hObject, 'merstruct', merstruct);

% merhandles is a struct array, each element has .side and .label (strings)
% for indexing, and .h (plotted object handle)
merhandles = struct('color_list', cat(1, options.prefs.mer.tract_info.color),...
    'traj', struct('side', {}, 'label', {}, 'h', []),...
    'markers', struct('side', {}, 'label', {}, 'depth', [], 'h', [], 'tag', []));
setappdata(hObject, 'merhandles', merhandles);

ea_gui_generate(handles);  % Generate UI elements for each tract.
ea_resultfig_clear(handles);  % Delete markers and mer trajectories.
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
[Depth,Thresh,Tracts] = ea_readmer(fullfile(options.root,options.patientname,['Rec_',sidestr,'.mat']));
merlfp.(sidestr).Depth = Depth;
merlfp.(sidestr).Tracts = Tracts;
merlfp.(sidestr).Thresh = Thresh;

save(fullfile(options.root,options.patientname,'ea_recordings.mat'),'merlfp')
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
ea_data_clear(handles);  % Clears trajectories, depths, etc in merstruct.
ea_resultfig_clear(handles);  % Delete markers and mer trajectories.
ea_mercontrol_updateall(handles);  % Set GUI from data.


% --- Executes during object creation, after setting all properties.
function clearall_CreateFcn(hObject, eventdata, handles)
% hObject    handle to clearall (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
background = imread(fullfile(ea_getearoot, 'icons/delete.png'));
set(hObject,'CData',background);


% --- Executes on button press in setdefault.
function setdefault_Callback(hObject, eventdata, handles)
% hObject    handle to setdefault (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ea_resultfig_clear(handles);  % Delete markers and mer trajectories visualizations.
merstruct = getappdata(handles.mercontrolfig, 'merstruct');
merstruct.setDataToDefaults();
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
background = imread(fullfile(ea_getearoot, 'icons/checkmark.png'));
set(hObject,'CData',background)


% --- Executes on button press in undomarker.
function undomarker_Callback(hObject, eventdata, handles)
% hObject    handle to undomarker (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
merstruct = getappdata(handles.mercontrolfig, 'merstruct');
merstruct.undoMarker();
ea_resultfig_updatemarkers(handles);  % Sync handles to match markers.
ea_mercontrol_updatemarkers(handles); % Enable/Disable marker gui elements.


% --- Executes during object creation, after setting all properties.
function undomarker_CreateFcn(hObject, eventdata, handles)
% hObject    handle to undomarker (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
background = imread(fullfile(ea_getearoot, 'icons/undo.png'));
set(hObject,'CData',background,'Value',0)


% --- Executes on button press in redomarker.
function redomarker_Callback(hObject, eventdata, handles)
% hObject    handle to redomarker (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
merstruct = getappdata(handles.mercontrolfig, 'merstruct');
merstruct.redoMarker();
ea_resultfig_updatemarkers(handles);  % Sync handles to match markers.
ea_mercontrol_updatemarkers(handles); % Enable/Disable marker gui elements.


% --- Executes during object creation, after setting all properties.
function redomarker_CreateFcn(hObject, eventdata, handles)
% hObject    handle to redomarker (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
background = imread(fullfile(ea_getearoot, 'icons/redo.png'));
set(hObject,'CData',background,'Value',0)


% --- Executes on button press in togglemarkertags.
function togglemarkertags_Callback(hObject, eventdata, handles)
% hObject    handle to togglemarkertags (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglemarkertags
merstruct = getappdata(handles.mercontrolfig, 'merstruct');
if get(hObject,'Value')
    merstruct.Config.vis.tag_visible = 'on';
else
    merstruct.Config.vis.tag_visible = 'off';
end
ea_resultfig_updatemarkers(handles);


% --- Executes during object creation, after setting all properties.
function togglemarkertags_CreateFcn(hObject, eventdata, handles)
% hObject    handle to togglemarkertags (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
background = imread(fullfile(ea_getearoot, 'icons/text.png'));
set(hObject,'CData',background,'Value',0)


% --- Executes on button press in loadstate.
function loadstate_Callback(hObject, eventdata, handles)
% hObject    handle to loadstate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
merstruct = getappdata(handles.mercontrolfig, 'merstruct');
merstruct.load();
ea_resultfig_update(handles);
ea_mercontrol_updateall(handles);


% --- Executes during object creation, after setting all properties.
function loadstate_CreateFcn(hObject, eventdata, handles)
% hObject    handle to loadstate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
background = imread(fullfile(ea_getearoot, 'icons/import.png'));
set(hObject,'CData',background,'Value',0)


% --- Executes on button press in savestate.
function savestate_Callback(hObject, eventdata, handles)
% hObject    handle to savestate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
merstruct = getappdata(handles.mercontrolfig, 'merstruct');
merstruct.save();


% --- Executes during object creation, after setting all properties.
function savestate_CreateFcn(hObject, eventdata, handles)
% hObject    handle to savestate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
background = imread(fullfile(ea_getearoot, 'icons/export.png'));
set(hObject,'CData',background,'Value',0)


% --- Executes on button press in pushbutton_clear_checks.
function pushbutton_clear_checks_Callback(hObject, eventdata, handles)
merstruct = getappdata(handles.mercontrolfig, 'merstruct');
merstruct.setTogglesToDefaults();
ea_mercontrol_updatetoggles(handles);


% --- Executes on button press in nexframeadj.
function nexframeadj_Callback(hObject, eventdata, handles)
% hObject    handle to nexframeadj (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Modal window for nexdrive rotation.
ea_nexframegui(handles.mercontrolfig);













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

merstruct = getappdata(handles.mercontrolfig, 'merstruct');
% Generate UI elements

width = handles.mainuipanel.Position(3)-5;

for side_str = {'left', 'right'}
    bSide = strcmpi({merstruct.Config.MERTrajectory.side}, side_str);
    side_tract_info = merstruct.Config.MERTrajectory(bSide);
    
    % In mainuipanel...
    
    % ...Add options to popupimplantedtract_left _right
    uq_labels = unique({side_tract_info.label});
    popup_string = cat(2, 'Set implanted...', uq_labels);
    set(handles.(['popupimplantedtract_' side_str{:}]),...
        'String', popup_string);
    
    % ...Add togglebuttons next to texttoggle_left and _right (with callbacks)
    anchor_pos = get(handles.(['texttoggle_' side_str{:}]), 'Position');
    w_per = (width - anchor_pos(3) - anchor_pos(1)) / length(side_tract_info);
    
    for tract_ix = 1:length(side_tract_info)
        pos_str = side_tract_info(tract_ix).label;
        new_pos = [anchor_pos(1) + anchor_pos(3) + (tract_ix-1)*w_per...
            anchor_pos(2)-5 w_per 26];
        uicontrol(handles.mainuipanel,...
            'Style', 'togglebutton',...
            'Tag', ['toggle' pos_str '_' side_str{:}],...
            'String', pos_str,...
            'Position', new_pos,...
            'Units', 'pixels',...
            'BackgroundColor', [1 1 1],...
            'Value', 1,...
            'Callback', @ea_updatetogglestate);
    end
    
    % ...Add edit boxes next to textpos_left and _right (with callbacks)
    anchor_pos = get(handles.(['textpos_' side_str{:}]), 'Position');
    w_per = (width - anchor_pos(3) - anchor_pos(1)) / length(side_tract_info);
    for tract_ix = 1:length(side_tract_info)
        pos_str = side_tract_info(tract_ix).label;
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
    w_per = (width - anchor_pos(3) - anchor_pos(1)) / length(side_tract_info);
    for tract_ix = 1:length(side_tract_info)
        pos_str = side_tract_info(tract_ix).label;
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
[side_strs, sid, track] = ea_detsidestr(hObject.Tag);
% Update the value
handles = guidata(hObject);
merstruct = getappdata(handles.mercontrolfig, 'merstruct');
bSide = strcmpi({merstruct.Toggles.keycontrol.side}, side_strs{sid});
bLabel = strcmpi({merstruct.Toggles.keycontrol.label}, track);
merstruct.Toggles.keycontrol(bSide & bLabel).value = get(hObject,'Value') == get(hObject,'Max');


function ea_data_clear(handles, varargin)
% ea_clear_data(handles)
% Resets trajectories, depths, etc.
merstruct = getappdata(handles.mercontrolfig, 'merstruct');
merstruct.setDataToDefaults();


function ea_resultfig_clear(handles)
% ea_resultfig_clear(handles)
% Delete markers and plot empty 3dplot trajs.
ea_resultfig_clear_traj(handles);
ea_resultfig_clear_markers(handles);


function ea_resultfig_clear_traj(handles)
resultfig = getappdata(handles.mercontrolfig, 'resultfig');
merstruct = getappdata(handles.mercontrolfig, 'merstruct');
merhandles = getappdata(handles.mercontrolfig, 'merhandles');

if isvalid(resultfig)
    % Get a list of tags for all the plotted items in resultfig.
    curr_fig_items = resultfig.CurrentAxes.Children;
    tag_array = arrayfun(@(x)x.Tag, curr_fig_items, 'UniformOutput', false)';
    
    set(0, 'CurrentFigure', resultfig); hold on;
    uqlabels = unique({merstruct.MERTrajectories.label}, 'stable');
    for traj_ix = 1:length(merstruct.MERTrajectories)
        traj = merstruct.MERTrajectories(traj_ix);
        bH = strcmpi({merhandles.traj.side}, traj.side) & strcmpi({merhandles.traj.label}, traj.label);
        delete([merhandles.traj(bH).h]);
        merhandles.traj(bH) = [];
        resid_items = curr_fig_items(strcmpi(tag_array, [traj.label '_' traj.side]));
        if ~isempty(resid_items)
            delete(resid_items);
        end
        % Create the plot handle.
        % Because we aren't plotting anything, the object handle will remain
        % empty. color and tag are ignored.
        this_color = merhandles.color_list(strcmpi(uqlabels, traj.label), :);
        h = plot3([], [], [], 'color', this_color,...
            'linew', 5, 'tag', [traj.label '_' traj.side]);
        merhandles.traj(end+1) = struct('side', traj.side, 'label', traj.label, 'h', h);
    end
    setappdata(handles.mercontrolfig, 'merhandles', merhandles);
end


function ea_resultfig_clear_markers(handles)
merhandles = getappdata(handles.mercontrolfig, 'merhandles');
for marker_ix = 1:length(merhandles.markers)
    try delete(merhandles.markers(marker_ix).h); end
    try delete(merhandles.markers(marker_ix).tag); end
    merhandles.markers(marker_ix).h = [];
    if isfield(merhandles.markers(marker_ix), 'tag')
        merhandles.markers(marker_ix).tag = [];
    end
end
setappdata(handles.mercontrolfig, 'merhandles', merhandles)


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


function ea_mercontrol_updatetoggles(handles)
merstruct = getappdata(handles.mercontrolfig, 'merstruct');
ui_tags = {handles.mainuipanel.Children.Tag};
kcg_tags = {handles.keyboardcontrolgroup.Children.Tag};
for ts_ix = 1:length(merstruct.Toggles.togglestates)
    ts = merstruct.Toggles.togglestates(ts_ix);
    set(handles.mainuipanel.Children(strcmpi(ui_tags,...
        ['toggle', ts.label, '_', ts.side])), 'Value', ts.value);
end
for kc_ix = 1:length(merstruct.Toggles.keycontrol)
    kc = merstruct.Toggles.keycontrol(kc_ix);
    set(handles.keyboardcontrolgroup.Children(strcmpi(kcg_tags,...
        ['keycheck', kc.label, '_', kc.side])), 'Value', kc.value);
end



function editimplanteddepth_helper(hObject, eventdata, handles)
depth = str2double(get(hObject, 'String'));
str_parts = strsplit(eventdata.Source.Tag, '_');
[sidestr, side_ix, ~] = ea_detsidestr(str_parts{2});
merstruct = getappdata(handles.mercontrolfig, 'merstruct');
merstruct.updateDBSDepth(sidestr{side_ix}, depth);
ea_mercontrol_updateimplanted(handles, sidestr{side_ix});
ea_resultfig_updatetrajectories(handles, sidestr{side_ix});


function popupimplantedtract_helper(hObject, eventdata, handles) %#ok<INUSL>
str_parts = strsplit(hObject.Tag, '_');
[sidestr, side_ix, ~] = ea_detsidestr(str_parts{2});
merstruct = getappdata(handles.mercontrolfig, 'merstruct');
merstruct.updateDBSImplantTrack(sidestr{side_ix}, hObject.String{hObject.Value});
ea_resultfig_updatetrajectories(handles, sidestr{side_ix});


function ea_updatepostext(hObject, eventData) %#ok<INUSD>
% Called when one of the position text boxes is edited.
handles = guidata(hObject);
[side_strs, sid, track] = ea_detsidestr(hObject.Tag);
merstruct = getappdata(handles.mercontrolfig, 'merstruct');
merstruct.updateTrajDepth(side_strs{sid}, track, str2double(get(hObject, 'String')));
ea_resultfig_updatetrajectories(handles);


function ea_updatetogglestate(hObject, eventdata) %#ok<INUSD>
% Called when one of the tract toggle buttons is pressed.
[side_strs, sid, track] = ea_detsidestr(hObject.Tag);
handles = guidata(hObject);
merstruct = getappdata(handles.mercontrolfig, 'merstruct');
bToggle = strcmpi({merstruct.Toggles.togglestates.side}, side_strs{sid}) ...
    & strcmpi({merstruct.Toggles.togglestates.label}, track);
merstruct.Toggles.togglestates(bToggle).value = hObject.Value;
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
