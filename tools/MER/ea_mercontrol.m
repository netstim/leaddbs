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

% Last Modified by GUIDE v2.5 19-Jun-2017 16:39:43

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


% --- Executes just before ea_mercontrol is made visible.
function ea_mercontrol_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)
clc
set(hObject,'Name','MER Control','CloseRequestFcn',@closefunciton,'KeyPressFcn',@ea_keypress)
% Choose default command line output for ea_mercontrol
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ea_mercontrol wait for user response (see UIRESUME)
% uiwait(handles.mercontrolfig);

resultfig=varargin{1};
options=varargin{2};
setappdata(hObject,'resultfig',resultfig);
setappdata(hObject,'options',options);

merstruct = options.prefs.mer;  %.length [24 mm]; .offset [2 mm]; .tract_info
setappdata(resultfig,'merstruct',merstruct)

% Generate UI elements in mainuipanel for each merstruct.tract_info
width = handles.mainuipanel.InnerPosition(3);
popup_string = cat(2, 'Set implanted...', {merstruct.tract_info.label});
for side_str = {'left', 'right'}
    % Add options to popupimplantedtract_left _right
    set(handles.(['popupimplantedtract_' side_str{:}]),...
        'String', popup_string);
    
    % Add togglebuttons next to texttoggle_left and _right (with callbacks)
    anchor_pos = get(handles.(['texttoggle_' side_str{:}]), 'Position');
    w_per = (width - anchor_pos(3) - anchor_pos(1)) / length(merstruct.tract_info);
    for tract_ix = 1:length(merstruct.tract_info)
        pos_str = merstruct.tract_info(tract_ix).label;
        new_pos = anchor_pos;
        new_pos(1) = anchor_pos(1) + anchor_pos(3) + (tract_ix-1)*w_per;
        new_pos(3) = w_per;
        htoggle = uicontrol(handles.mainuipanel,...
            'Style', 'togglebutton',...
            'Tag', ['toggle' pos_str '_' side_str{:}],...
            'Value', 1,...
            'String', pos_str,...
            'Position', new_pos,...
            'Units', 'pixels',...
            'BackgroundColor', [1 1 1],...
            'Value', 0,...
            'Callback', @ea_updatetogglestate);
    end
    
    % Add edit boxes next to textpos_left and _right (with callbacks)
    anchor_pos = get(handles.(['textpos_' side_str{:}]), 'Position');
    w_per = (width - anchor_pos(3) - anchor_pos(1)) / length(merstruct.tract_info);
    for tract_ix = 1:length(merstruct.tract_info)
        pos_str = merstruct.tract_info(tract_ix).label;
        new_pos = anchor_pos;
        new_pos(1) = anchor_pos(1) + anchor_pos(3) + (tract_ix-1)*w_per;
        new_pos(3) = w_per;
        hedit = uicontrol(handles.mainuipanel,...
            'Style', 'edit',...
            'Tag', ['pos' pos_str '_' side_str{:}],...
            'String', '0',...
            'Position', new_pos,...
            'Units', 'pixels',...
            'BackgroundColor', [1 1 1],...
            'Callback', @ea_updatepostext);
    end
    
    % In keyboardcontrolgroup, add radiobuttons next to keytext_left and _right inside
    anchor_pos = get(handles.(['textkey_' side_str{:}]), 'Position');
    w_per = (width - anchor_pos(3) - anchor_pos(1)) / length(merstruct.tract_info);
    for tract_ix = 1:length(merstruct.tract_info)
        pos_str = merstruct.tract_info(tract_ix).label;
        new_pos = anchor_pos;
        new_pos(1) = anchor_pos(1) + anchor_pos(3) + (tract_ix-1)*w_per;
        new_pos(3) = w_per;
        hcheck = uicontrol(handles.keyboardcontrolgroup,...
            'Style', 'checkbox',...
            'Tag', ['keycheck' pos_str '_' side_str{:}],...
            'String', pos_str,...
            'Position', new_pos,...
            'Units', 'pixels',...
            'BackgroundColor', [1 1 1],...
            'Callback', @ea_keycheck);
    end
end
clearkeyboardcontrolgroup(handles);
clearmertrajectories(handles,resultfig,options);

[~,trajectory,markers]=ea_load_reconstruction(options);
coords_mm=ea_resolvecoords(markers,options);

% Set mermarkers
mermarkers = getappdata(resultfig,'mermarkers');
if isempty(mermarkers)
    mermarkers = struct('side',{},'tract',{},'depth',{},'markertype',{},...
        'session',{},'dat',{},'tag',{},'handle',{},'notes',{});
    setappdata(resultfig,'mermarkers',mermarkers)
end
% Set default implanted tract
set(handles.popupimplantedtract_left,'Value',options.prefs.mer.defaulttract+1)
set(handles.popupimplantedtract_right,'Value',options.prefs.mer.defaulttract+1)
set(handles.editimplanteddepth_left,'String','0')
set(handles.editimplanteddepth_right,'String','0')

merstruct.defaultmer = ea_coords2mer(coords_mm,trajectory,resultfig,options);
clear coords_mm trajectory markers

set(0,'CurrentFigure',resultfig); hold on;
for side = options.sides
    % contact_spacing = getfield(getappdata(resultfig,'elspec'),'contact_spacing');
    for pos_ix = 1:length(merstruct.tract_info)
        pos_str = merstruct.tract_info(pos_ix).label;
        merstruct.defaultmer.(pos_str).trajectory{side} = ...
            ea_getmertrajectory(merstruct.defaultmer.(pos_str).coords_mm{side}, 0, merstruct.length, 50);
        traj = merstruct.defaultmer.(pos_str).trajectory{side};
        merhandles.(pos_str){side} = plot3(traj(:,1), traj(:,2), traj(:,3),...
            'color',merstruct.tract_info(pos_ix).color, 'linew',5);
    end
end

merstruct.currentmer=merstruct.defaultmer;
setappdata(resultfig,'merstruct',merstruct)
setappdata(resultfig,'merhandles',merhandles)
ea_readtogglestates(handles);


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
% hObject    handle to popupimplantedtract_left (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupimplantedtract_left contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupimplantedtract_left
ea_set_implanted_helper(handles, 2);


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
% hObject    handle to popupimplantedtract_right (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupimplantedtract_right contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupimplantedtract_right
ea_set_implanted_helper(handles, 1);


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
% hObject    handle to editimplanteddepth_left (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editimplanteddepth_left as text
%        str2double(get(hObject,'String')) returns contents of editimplanteddepth_left as a double
side = 2; %left
dist = str2double(get(hObject,'String'));
resultfig=getappdata(handles.mercontrolfig,'resultfig');
options=getappdata(handles.mercontrolfig,'options');
merhandles = getappdata(resultfig,'merhandles');
refreshresultfig(handles,resultfig,options,side)

set(handles.editimplanteddepth_left,'String',num2str(dist));
% ea_getsettogglestates(handles,side,dist);


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
% hObject    handle to editimplanteddepth_right (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editimplanteddepth_right as text
%        str2double(get(hObject,'String')) returns contents of editimplanteddepth_right as a double
side = 1; %right
dist = str2double(get(hObject,'String'));
resultfig=getappdata(handles.mercontrolfig,'resultfig');
options=getappdata(handles.mercontrolfig,'options');
refreshresultfig(handles,resultfig,options,side)

set(handles.editimplanteddepth_right,'String',num2str(dist));
% ea_getsettogglestates(handles,side,dist);


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


function ea_keycheck(hObject, eventdata)
handles = guidata(hObject);
resultfig = getappdata(handles.mercontrolfig, 'resultfig');
merstruct = getappdata(resultfig, 'merstruct');
mertoggles = getappdata(resultfig, 'mertoggles');

% Determine which button called this
str_parts = split(hObject.Tag(9:end), '_');
tract_ix = find(strcmpi({merstruct.tract_info.label}, str_parts{1}));
side_ix = 1;
if strcmpi(str_parts{2}, 'left')
    side_ix = 2;
end

% Update the value
mertoggles.keycontrol(side_ix, tract_ix) = get(hObject,'Value') == get(hObject,'Max');
setappdata(resultfig, 'mertoggles', mertoggles);


% --- Executes when selected object is changed in keyboardcontrolgroup.
function keyboardcontrolgroup_SelectionChangedFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in keyboardcontrolgroup 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% mertoggles=getappdata(getappdata(gcf,'resultfig'),'mertoggles');
keymer = get(eventdata.NewValue,'tag');
% switch get(eventdata.NewValue,'tag')
%     case 'keycentral_left'
%            
%     case 'keyposterior_left'
%         
%     case 'keymedial_left'
%         
%     case 'keycentral_right'
% 
%     case 'keyposterior_right'
%         
%     case 'keymedial_right'
%         
% end
setappdata(getappdata(gcf,'resultfig'),'keymer',keymer);


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

% Hint: get(hObject,'Value') returns toggle state of clearall
resultfig = getappdata(handles.mercontrolfig,'resultfig');
options = getappdata(handles.mercontrolfig,'options');
clearmertrajectories(handles,resultfig,options)
clearmermarkers(handles,resultfig,options)
set(hObject,'Value',0)


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
set(hObject,'CData',background,'Value',0)


% --- Executes on button press in setdefault.
function setdefault_Callback(hObject, eventdata, handles)
% hObject    handle to setdefault (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of setdefault
resultfig = getappdata(handles.mercontrolfig,'resultfig');
options = getappdata(handles.mercontrolfig,'options');
merstruct = getappdata(resultfig,'merstruct');
setdefaultmer(handles,resultfig,merstruct,options)
ea_getsettogglestates(hObject, eventdata, handles);
set(hObject,'Value',0)


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
set(hObject,'CData',background,'Value',0)


% --- Executes on button press in undomarker.
function undomarker_Callback(hObject, eventdata, handles)
% hObject    handle to undomarker (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of undomarker
resultfig = getappdata(handles.mercontrolfig,'resultfig');
mermarkers = getappdata(resultfig,'mermarkers');
markerstring = getappdata(resultfig,'markerstring');
n = length(mermarkers);
if n~=0
    tmpmer = mermarkers(n);
    mermarkers(n)=[];
    set(tmpmer.handle,'Visible','off')
    set(tmpmer.tag.handle,'Visible','off')
    setappdata(resultfig,'tmpmer',tmpmer);
    
    if strcmpi(tmpmer.side,'right')
        markerstring.right(end) = [];
    elseif strcmpi(tmpmer.side,'left')
        markerstring.left(end) = [];
    end
    
    setappdata(resultfig,'mermarkers',mermarkers);
    setappdata(resultfig,'markerstring',markerstring)
    set(handles.popupmermarkers_right,'Visible','on','String',markerstring.right,'Value',1)
    set(handles.popupmermarkers_left,'Visible','on','String',markerstring.left,'Value',1)

end
set(hObject,'Value',0)


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

% Hint: get(hObject,'Value') returns toggle state of redomarker
resultfig = getappdata(handles.mercontrolfig,'resultfig');
options = getappdata(handles.mercontrolfig,'options');
mermarkers = getappdata(resultfig,'mermarkers');
markerstring = getappdata(resultfig,'markerstring');
tmpmer = getappdata(resultfig,'tmpmer');
n = length(mermarkers);
set(0,'CurrentFigure',resultfig);
if ~isempty(tmpmer)
    mermarkers(n+1)=tmpmer;
    % display marker
    tmp = tmpmer.handle;
    mermarkers(n+1).handle  = surf(tmp.XData,tmp.YData,tmp.ZData,'FaceLighting','gouraud',...
        'FaceColor',tmp.FaceColor,'EdgeColor',tmp.EdgeColor,'FaceAlpha',tmp.FaceAlpha,'Visible','on');
    clear tmp
    
    % display tag
    tmp = tmpmer.tag.handle;
    mermarkers(n+1).tag.handle = text(tmp.Position(1),tmp.Position(2),tmp.Position(3),...
        tmp.String,'Color',tmp.Color,'HorizontalAlignment',tmp.HorizontalAlignment,'Visible','on');
    clear tmp
    
    keymer = tmpmer.name;
    if strcmp(keymer(strfind(keymer,'_')+1:end),'right')
        markerstring.right{end+1} = sprintf('%0.0f. %s',n+1,tmpmer.tag.string);
    elseif strcmp(keymer(strfind(keymer,'_')+1:end),'left')
        markerstring.left{end+1} = sprintf('%0.0f. %s',n+1,tmpmer.tag.string);
    end
    
    tmp = []; %rmfield(getappdata(resultfig),'tmpmer');
    setappdata(resultfig,'tmpmer',tmp);
    setappdata(resultfig,'mermarkers',mermarkers);
    setappdata(resultfig,'markerstring',markerstring)
    set(handles.popupmermarkers_right,'Visible','on','String',markerstring.right,'Value',1)
    set(handles.popupmermarkers_left,'Visible','on','String',markerstring.left,'Value',1)

end
set(hObject,'Value',0)


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
options = getappdata(handles.mercontrolfig,'options');
mermarkers = getappdata(resultfig,'mermarkers');
if get(hObject,'Value')
    for i = 1:length(mermarkers)
        try set(mermarkers(i).tag.handle,'Visible','on'); end
    end
    % set(handles.togglemarkertags,'Value',1)
else
    for i = 1:length(mermarkers)
        try set(mermarkers(i).tag.handle,'Visible','off'); end
    end
    % set(handles.togglemarkertags,'Value',0)
end


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

% Hint: get(hObject,'Value') returns toggle state of importmarkers
resultfig = getappdata(handles.mercontrolfig,'resultfig');
options = getappdata(handles.mercontrolfig,'options');
merhandles = getappdata(resultfig,'merhandles');
side=1;
if isvalid(merhandles.central{side})
try 
    load([options.uipatdirs{1},filesep,'ea_mermarkers.mat'],'mermarkers','mertoggles');
catch
    [files] = fdir(options.uipatdirs{1},'.mat');
    if ~isempty(cell2mat(strfind({files(:).name},'ea_mermarkers')))
        %[s,b] = listdlg('PromptString','Select a file:','SelectionMode','single',...
        %   'ListString',{files(~cellfun(@isempty,strfind({files(:).name},'ea_mermarkers'))).name})
    end
    ea_error(['could not find ea_mermarkers.mat for ', options.patientname])
end
%assignin(workspace,'mermarkers',mermarkers)
ea_getsettogglestates(hObject, [], handles, mertoggles);
ea_updatemermarkers(mermarkers,handles,resultfig,options)
else
    fprintf(2,'Error using importmarkers_Callback: \nNo MER trajectories available\n\n')
end
set(hObject,'Value',0)


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

% Hint: get(hObject,'Value') returns toggle state of exportmarkers
resultfig = getappdata(handles.mercontrolfig,'resultfig');
options = getappdata(handles.mercontrolfig,'options');
mermarkers = getappdata(resultfig,'mermarkers');
mertoggles = ea_getsettogglestates(hObject, [], handles);
if exist([options.uipatdirs{1},filesep,'ea_mermarkers.mat'],'file')
    overwrite = ea_questdlg({['ea_mermarkers.mat found in ' options.patientname ' directory.'];...
        'Would you like to overwrite this file?'},options.patientname);
end
if ~exist([options.uipatdirs{1},filesep,'ea_mermarkers.mat'],'file') || strcmp(overwrite,'Yes')
    disp(['Saving: ' options.uipatdirs{1},filesep,'ea_mermarkers.mat'])
    save([options.uipatdirs{1},filesep,'ea_mermarkers.mat'],'mermarkers','mertoggles')
    disp('DONE')
end    
set(hObject,'Value',0)


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




% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
function setdefaultmer(handles,resultfig,merstruct,options,varargin)
if isempty(varargin) && isequal(options.sides,1:2)
    sidestr='both';
elseif ~isempty(varargin)
    sidestr=varargin{1};
end
if strcmpi(sidestr, 'both')
    clear_sides = 1:2;
elseif strcmpi(sidestr, 'right')
    clear_sides = 1;
else
    clear_sides = 2;  % left
end
if length(varargin) > 1
    implanteddepth=varargin{2};
end
clearmertrajectories(handles,resultfig,options,sidestr)
merhandles=getappdata(resultfig,'merhandles');
set(0,'CurrentFigure',resultfig)

for side = clear_sides
    for pos_ix = 1:length(merstruct.tract_labels)
        pos_str = merstruct.tract_labels{pos_ix};
        traj = merstruct.defaultmer.(pos_str).trajectory{side};
        merhandles.(pos_str){side} = plot3(traj(:,1), traj(:,2), traj(:,3),...
            'color',merstruct.colormap(pos_ix,:), 'linew',5);
    end
end

setappdata(resultfig,'merhandles',merhandles)


function clearmermarkers(handles,resultfig,options)
% clear mermarkers and reset structure
mermarkers = getappdata(resultfig,'mermarkers');
for i = 1:length(mermarkers)
try delete(mermarkers(i).handle); end
try delete(mermarkers(i).tag.handle); end
end
mermarkers = struct('side',{},'tract',{},'depth',{},'markertype',{},'session',{},'dat',{},'tag',{},'handle',{},'notes',{});
setappdata(resultfig,'mermarkers',mermarkers)
set(handles.popupimplantedtract_left,'Value',options.prefs.mer.defaulttract+1);
set(handles.popupimplantedtract_right,'Value',options.prefs.mer.defaulttract+1);
set(handles.editimplanteddepth_left,'String','0');
set(handles.editimplanteddepth_right,'String','0');
set(handles.popupmermarkers_left,'Visible','off','String','','Value',1)
set(handles.popupmermarkers_right,'Visible','off','String','','Value',1)


function clearmertrajectories(handles,resultfig,options,sidestr)

merhandles = getappdata(resultfig,'merhandles');
if isempty(merhandles)
    return
end
if ~exist('sidestr','var') && isequal(options.sides,1:2)
    clear_sides = 1:2;
% Legacy
elseif isnumeric(sidestr)
    clear_sides = sidestr;
else
    if strcmpi(sidestr, 'both')
        clear_sides = 1:2;
    elseif strcmpi(sidestr, 'right')
        clear_sides = 1;
    elseif strcmpi(sidestr, 'left')
        clear_sides = 2;
    end
end

merstruct = getappdata(resultfig,'merstruct');
ui_tags = {handles.mainuipanel.Children.Tag};
for side = clear_sides
    if side == 1
        sidestr = 'right';
    elseif side == 2
        sidestr = 'left';
    end
    for pos_ix = 1:length(merstruct.tract_info)
        pos_str = merstruct.tract_info(pos_ix).label;
        try; delete(merhandles.(pos_str){side}); end
        edit_handle = handles.mainuipanel.Children(strcmpi(ui_tags, ['pos' pos_str '_' sidestr]));
        set(edit_handle, 'String', '0');
    end
end
set(handles.editimplanteddepth_left,'String','0');
set(handles.editimplanteddepth_right,'String','0');
setappdata(resultfig,'merhandles',merhandles)


function closefunciton(hObject,event)
try
resultfig = getappdata(hObject,'resultfig');
options = getappdata(hObject,'options');
handles = getappdata(hObject,'UsedByGUIData_m');
clearmertrajectories(handles,resultfig,options)
clearmermarkers(handles,resultfig,options)
end
delete(hObject)

function ea_updatemermarkers(mermarkers,handles,resultfig,options)
clearmermarkers(handles,resultfig,options)
set(0,'CurrentFigure',resultfig)
markerstring.right = {'none selected...'};
markerstring.left = {'none selected...'};
for n = 1:length(mermarkers)
    tmp = mermarkers(n).handle;
    mermarkers(n).handle = surf(tmp.XData,tmp.YData,tmp.ZData,...
        'FaceColor',tmp.FaceColor,'EdgeColor',tmp.EdgeColor,'FaceAlpha',tmp.FaceAlpha); 
    clear tmp
    tmp = mermarkers(n).tag.handle;
    mermarkers(n).tag.handle = text(tmp.Position(1),tmp.Position(2),tmp.Position(3),...
        tmp.String,'Color',tmp.Color,'HorizontalAlignment',tmp.HorizontalAlignment,'Visible','off');

    markerstring.(mermarkers(n).side){end+1} = sprintf('%0.0f. %s',n,mermarkers(n).tag.string);    
end

% set implanted tract data
try
    Lidx = find(~cellfun(@isempty,strfind({mermarkers.side},'left'))==1);
    set(handles.popupimplantedtract_left,'Value',...
        find(~cellfun(@isempty,strfind(handles.popupimplantedtract_left.String,mermarkers(Lidx(1)).dat.implantedtract))==1))
    set(handles.editimplanteddepth_left,'String',num2str(mermarkers(Lidx(1)).dat.leaddepth))
    ea_getsettogglestates(hObject, [], handles, 2, mermarkers(Lidx(1)).dat.leaddepth);
end

try
    Ridx = ~cellfun(@isempty,strfind({mermarkers.side},'right'));
    set(handles.popupimplantedtract_right,'Value',...
        find(~cellfun(@isempty,strfind(handles.popupimplantedtract_left.String,mermarkers(Ridx(1)).dat.implantedtract))==1))
    set(handles.editimplanteddepth_left,'String',num2str(mermarkers(Ridx(1)).dat.leaddepth))
    ea_getsettogglestates(hObject, [], handles, 1, mermarkers(Ridx(1)).dat.leaddepth);
end

setappdata(resultfig,'mermarkers',mermarkers)
setappdata(resultfig,'markerstring',markerstring)
%set(handles.popupmermarkers_right,'Visible','on','String',markerstring.right,'Value',1)
%set(handles.popupmermarkers_left,'Visible','on','String',markerstring.left,'Value',1)
set(handles.togglemarkertags,'Visible','on','Value',0)


function ea_updatepostext(hObject, eventData)
% Called when one of the position text boxes is edited.
handles = guidata(hObject);
resultfig = getappdata(handles.mercontrolfig, 'resultfig');
merstruct = getappdata(resultfig, 'merstruct');

% Determine which text element called this
str_parts = split(hObject.Tag(4:end), '_');
% tract_ix = find(strcmpi({merstruct.tract_info.label}, str_parts{1}));
side_ix = 1;
if strcmpi(str_parts{2}, 'left')
    side_ix = 2;
end

input = str2double(get(hObject,'String'));
implant_depth = str2double(get(handles.(['editimplanteddepth_' str_parts{2}]),'String'));
dist = input - implant_depth;
trajectory = ea_getmertrajectory(merstruct.currentmer.(str_parts{1}).coords_mm{side_ix},...
    dist, merstruct.length, 50);
ea_updatemertrajectory(handles,...
    trajectory, 0, hObject.Tag);



function mertoggles = ea_updatetogglestate(hObject, eventdata)
% Called when one of the tract toggle buttons is pressed.
handles = guidata(hObject);
resultfig = getappdata(handles.mercontrolfig, 'resultfig');
merstruct = getappdata(resultfig, 'merstruct');
mertoggles = getappdata(resultfig, 'mertoggles');

% Determine which button called this
str_parts = split(hObject.Tag(7:end), '_');
tract_ix = find(strcmpi({merstruct.tract_info.label}, str_parts{1}));
side_ix = 1;
if strcmpi(str_parts{2}, 'left')
    side_ix = 2;
end

% Update the value in mertoggles
mertoggles.togglestates(side_ix, tract_ix) = hObject.Value;
setappdata(resultfig,'mertoggles',mertoggles);

% Update the visual presentation
merhandles = getappdata(resultfig,'merhandles');
if isvalid(merhandles.(str_parts{1}){side_ix})
    if hObject.Value
        set(merhandles.(str_parts{1}){side_ix}, 'Visible', 'on');
    else
        set(merhandles.(str_parts{1}){side_ix}, 'Visible', 'off');
    end
end


function mertoggles = ea_readtogglestates(handles)
% Scans the toggle button values and stores them in appdata
resultfig = getappdata(handles.mercontrolfig, 'resultfig');
merstruct = getappdata(resultfig, 'merstruct');
ui_tags = {handles.mainuipanel.Children.Tag};
side_strs = {'right', 'left'};
mertoggles.togglestates = nan(2, length(merstruct.tract_info));
for side_ix = 1:length(side_strs)
    for pos_ix = 1:length(merstruct.tract_info)
        search_str = ['toggle' merstruct.tract_info(pos_ix).label '_' side_strs{side_ix}];
        mertoggles.togglestates(side_ix, pos_ix) = handles.mainuipanel.Children(strcmp(ui_tags, search_str)).Value;
        
        % TODO: keycheck
    end
end
setappdata(resultfig,'mertoggles',mertoggles);


function mertoggles = ea_getsettogglestates(hObject, eventdata, handles, varargin)
% Example type A: 
%   varargin{1}=mertoggles;
%
% Example type B: ea_getsettogglestates(handles,1,mermarkers(Ridx(1)).dat.leaddepth);
%   varargin{1}=side;
%   varargin{2}=dist;

if size(varargin,2)==1
    mertoggles = varargin{1};
elseif size(varargin,2)==2
    side = varargin{1};
    dist = varargin{2};
end

resultfig=getappdata(handles.mercontrolfig,'resultfig');
merhandles = getappdata(resultfig,'merhandles');
options = getappdata(handles.mercontrolfig,'options');
merstruct = getappdata(resultfig, 'merstruct');

if ~exist('mertoggles', 'var')
    mertoggles.keycontrol = nan(2, length(merstruct.tract_info));
    mertoggles.togglestates = nan(2, length(merstruct.tract_info));
end

% reset states based on gui:
for side = 1:2
    if side == 1
        side_str = 'right';
    else
        side_str = 'left';
    end
    for pos_ix = 1:length(merstruct.tract_info)
        pos_str = merstruct.tract_info(pos_ix).label;
        this_substr = [pos_str '_' side_str];
        if isnan(mertoggles.keycontrol(side, pos_ix))
            mertoggles.keycontrol(side, pos_ix) = get(handles.(['key' this_substr]),'Value');
            mertoggles.togglestates(side, pos_ix) = get(handles.(['toggle' this_substr]),'Value');
        else
            set(handles.(['key' this_substr]),'Value',mertoggles.keycontrol(side, pos_ix));
            set(handles.(['toggle' this_substr]),'Value', mertoggles.togglestates(side, pos_ix));
        end
        % set togglestates 
        if get(handles.(['toggle' this_substr]),'Value') && isvalid(merhandles.(pos_str){side})
            set(merhandles.(pos_str){side},'Visible','on')
        elseif ~get(handles.(['toggle' this_substr]),'Value') && isvalid(merhandles.(pos_str){side})
            set(merhandles.(pos_str){side},'Visible','off')
        end
        % set pos handles
        if exist('dist', 'var')
            set(handles.(['pos' this_substr]),'String',num2str(dist));
        end
    end
end
setappdata(getappdata(handles.mercontrolfig,'resultfig'),'mertoggles',mertoggles); % also store toggle data in resultfig.
setappdata(resultfig,'merhandles',merhandles)


function merstruct = ea_updatetracts(merstruct,side,tracttag,resultfig,handles)
        
offset = getfield(getappdata(resultfig,'merstruct'),'offset');
% elstruct = 
% x-axis --> negative = medial, positive = lateral
% y-axis --> negative = posterior, positive = anterior
% z-axis --> negative = inferior, positive = superior

% defaultmer coords assumes central is implant site
coords_mm = merstruct.defaultmer.central.coords_mm{side};
trajectory = merstruct.defaultmer.central.trajectory{side};
switch tracttag 
    case 'central'
        merstruct.currentmer.central.coords_mm{side} = coords_mm;
        merstruct.currentmer.central.trajectory{side} = trajectory;
        merstruct.currentmer.anterior.trajectory{side} = [coords_mm(:,1),coords_mm(:,2)+offset,coords_mm(:,3)]; %2mm anterior
        merstruct.currentmer.anterior.trajectory{side} = [trajectory(:,1),trajectory(:,2)+offset,trajectory(:,3)]; %2mm anterior
        merstruct.currentmer.posterior.coords_mm{side} = [coords_mm(:,1),coords_mm(:,2)-offset,coords_mm(:,3)]; %2mm posterior
        merstruct.currentmer.posterior.trajectory{side} = [trajectory(:,1),trajectory(:,2)-offset,trajectory(:,3)]; %2mm posterior
        if side==1 %right
            merstruct.currentmer.lateral.coords_mm{side} = [coords_mm(:,1)+offset,coords_mm(:,2),coords_mm(:,3)]; %2mm lateral (right is positive)
            merstruct.currentmer.lateral.trajectory{side} = [trajectory(:,1)+offset,trajectory(:,2),trajectory(:,3)]; %2mm lateral (right is positive)
            merstruct.currentmer.medial.coords_mm{side} = [coords_mm(:,1)-offset,coords_mm(:,2),coords_mm(:,3)]; %2mm medial (right is positive)
            merstruct.currentmer.medial.trajectory{side} = [trajectory(:,1)-offset,trajectory(:,2),trajectory(:,3)]; %2mm medial (right is positive)
        elseif side==2 %left
            merstruct.currentmer.lateral.coords_mm{side} = [coords_mm(:,1)-offset,coords_mm(:,2),coords_mm(:,3)]; %2mm medial (left is negative)
            merstruct.currentmer.lateral.trajectory{side} = [trajectory(:,1)-offset,trajectory(:,2),trajectory(:,3)]; %2mm medial (left is negative)
            merstruct.currentmer.medial.coords_mm{side} = [coords_mm(:,1)+offset,coords_mm(:,2),coords_mm(:,3)]; %2mm medial (left is negative)
            merstruct.currentmer.medial.trajectory{side} = [trajectory(:,1)+offset,trajectory(:,2),trajectory(:,3)]; %2mm medial (left is negative)
        end
    case 'anterior'
        merstruct.currentmer.anterior.coords_mm{side} = coords_mm;
        merstruct.currentmer.anterior.trajectory{side} = trajectory;
        merstruct.currentmer.central.coords_mm{side} = [coords_mm(:,1),coords_mm(:,2)-offset,coords_mm(:,3)]; %2mm anterior
        merstruct.currentmer.central.trajectory{side} = [trajectory(:,1),trajectory(:,2)-offset,trajectory(:,3)]; %2mm anterior
        merstruct.currentmer.posterior.coords_mm{side} = [coords_mm(:,1),coords_mm(:,2)-offset*2,coords_mm(:,3)]; %4mm anterior
        merstruct.currentmer.posterior.trajectory{side} = [trajectory(:,1),trajectory(:,2)-offset*2,trajectory(:,3)]; %4mm anterior
         if side==1 %right
            merstruct.currentmer.lateral.coords_mm{side} =  [coords_mm(:,1)+offset,coords_mm(:,2)-offset,coords_mm(:,3)]; %2mm medial (right is positive)
            merstruct.currentmer.lateral.trajectory{side} = [trajectory(:,1)+offset,trajectory(:,2)-offset,trajectory(:,3)]; %2mm medial (right is positive)
            merstruct.currentmer.medial.coords_mm{side} = [coords_mm(:,1)-offset,coords_mm(:,2)-offset,coords_mm(:,3)]; %2mm medial (right is positive)
            merstruct.currentmer.medial.trajectory{side} = [trajectory(:,1)-offset,trajectory(:,2)-offset,trajectory(:,3)]; %2mm medial (right is positive)
        elseif side==2 %left
            merstruct.currentmer.lateral.coords_mm{side} = [coords_mm(:,1)-offset,coords_mm(:,2)-offset,coords_mm(:,3)]; %2mm medial (left is negative)
            merstruct.currentmer.lateral.trajectory{side} = [trajectory(:,1)-offset,trajectory(:,2)-offset,trajectory(:,3)]; %2mm medial (left is negative)
            merstruct.currentmer.medial.coords_mm{side} = [coords_mm(:,1)+offset,coords_mm(:,2)-offset,coords_mm(:,3)]; %2mm medial (left is negative)
            merstruct.currentmer.medial.trajectory{side} = [trajectory(:,1)+offset,trajectory(:,2)-offset,trajectory(:,3)]; %2mm medial (left is negative)
         end
    case 'posterior'
        merstruct.currentmer.posterior.coords_mm{side} = coords_mm;
        merstruct.currentmer.posterior.trajectory{side} = trajectory;
        merstruct.currentmer.central.coords_mm{side} = [coords_mm(:,1),coords_mm(:,2)+offset,coords_mm(:,3)]; %2mm anterior
        merstruct.currentmer.central.trajectory{side} = [trajectory(:,1),trajectory(:,2)+offset,trajectory(:,3)]; %2mm anterior
        merstruct.currentmer.anterior.coords_mm{side} = [coords_mm(:,1),coords_mm(:,2)+offset*2,coords_mm(:,3)]; %4mm anterior
        merstruct.currentmer.anterior.trajectory{side} = [trajectory(:,1),trajectory(:,2)+offset*2,trajectory(:,3)]; %4mm anterior
         if side==1 %right
            merstruct.currentmer.lateral.coords_mm{side} =  [coords_mm(:,1)+offset,coords_mm(:,2)+offset,coords_mm(:,3)]; %2mm medial (right is positive)
            merstruct.currentmer.lateral.trajectory{side} = [trajectory(:,1)+offset,trajectory(:,2)+offset,trajectory(:,3)]; %2mm medial (right is positive)
            merstruct.currentmer.medial.coords_mm{side} = [coords_mm(:,1)-offset,coords_mm(:,2)+offset,coords_mm(:,3)]; %2mm medial (right is positive)
            merstruct.currentmer.medial.trajectory{side} = [trajectory(:,1)-offset,trajectory(:,2)+offset,trajectory(:,3)]; %2mm medial (right is positive)
        elseif side==2 %left
            merstruct.currentmer.lateral.coords_mm{side} = [coords_mm(:,1)-offset,coords_mm(:,2)+offset,coords_mm(:,3)]; %2mm medial (left is negative)
            merstruct.currentmer.lateral.trajectory{side} = [trajectory(:,1)-offset,trajectory(:,2)+offset,trajectory(:,3)]; %2mm medial (left is negative)
            merstruct.currentmer.medial.coords_mm{side} = [coords_mm(:,1)+offset,coords_mm(:,2)+offset,coords_mm(:,3)]; %2mm medial (left is negative)
            merstruct.currentmer.medial.trajectory{side} = [trajectory(:,1)+offset,trajectory(:,2)+offset,trajectory(:,3)]; %2mm medial (left is negative)
         end
    case 'lateral'
        merstruct.currentmer.lateral.coords_mm{side} = coords_mm;
        merstruct.currentmer.lateral.trajectory{side} = trajectory;
        if side==1 %right
            merstruct.currentmer.central.coords_mm{side} = [coords_mm(:,1)-offset,coords_mm(:,2),coords_mm(:,3)];
            merstruct.currentmer.central.trajectory{side} = [trajectory(:,1)-offset,trajectory(:,2),trajectory(:,3)];
            merstruct.currentmer.anterior.coords_mm{side} = [coords_mm(:,1)-offset,coords_mm(:,2)+offset,coords_mm(:,3)];
            merstruct.currentmer.anterior.trajectory{side} = [trajectory(:,1)-offset,trajectory(:,2)+offset,trajectory(:,3)];
            merstruct.currentmer.posterior.coords_mm{side} = [coords_mm(:,1)-offset,coords_mm(:,2)-offset,coords_mm(:,3)];
            merstruct.currentmer.posterior.trajectory{side} = [trajectory(:,1)-offset,trajectory(:,2)-offset,trajectory(:,3)];
            merstruct.currentmer.medial.coords_mm{side} =  [coords_mm(:,1)-offset*2,coords_mm(:,2),coords_mm(:,3)]; %2mm medial (right is positive)
            merstruct.currentmer.medial.trajectory{side} = [trajectory(:,1)-offset*2,trajectory(:,2),trajectory(:,3)]; %2mm medial (right is positive)
        elseif side==2 %left
            merstruct.currentmer.central.coords_mm{side} = [coords_mm(:,1)+offset,coords_mm(:,2),coords_mm(:,3)];
            merstruct.currentmer.central.trajectory{side} = [trajectory(:,1)+offset,trajectory(:,2),trajectory(:,3)];
            merstruct.currentmer.anterior.coords_mm{side} = [coords_mm(:,1)+offset,coords_mm(:,2)+offset,coords_mm(:,3)];
            merstruct.currentmer.anterior.trajectory{side} = [trajectory(:,1)+offset,trajectory(:,2)+offset,trajectory(:,3)];
            merstruct.currentmer.posterior.coords_mm{side} = [coords_mm(:,1)+offset,coords_mm(:,2)-offset,coords_mm(:,3)];
            merstruct.currentmer.posterior.trajectory{side} = [trajectory(:,1)+offset,trajectory(:,2)-offset,trajectory(:,3)];
            merstruct.currentmer.medial.coords_mm{side} =  [coords_mm(:,1)+offset*2,coords_mm(:,2),coords_mm(:,3)]; %2mm medial (right is positive)
            merstruct.currentmer.medial.trajectory{side} = [trajectory(:,1)+offset*2,trajectory(:,2),trajectory(:,3)]; %2mm medial (right is positive)
        end
    case 'medial'
        merstruct.currentmer.medial.coords_mm{side} = coords_mm;
        merstruct.currentmer.medial.trajectory{side} = trajectory;
        if side==1 %right
            merstruct.currentmer.central.coords_mm{side} = [coords_mm(:,1)+offset,coords_mm(:,2),coords_mm(:,3)];
            merstruct.currentmer.central.trajectory{side} = [trajectory(:,1)+offset,trajectory(:,2),trajectory(:,3)];
            merstruct.currentmer.anterior.coords_mm{side} = [coords_mm(:,1)+offset,coords_mm(:,2)+offset,coords_mm(:,3)];
            merstruct.currentmer.anterior.trajectory{side} = [trajectory(:,1)+offset,trajectory(:,2)+offset,trajectory(:,3)];
            merstruct.currentmer.posterior.coords_mm{side} = [coords_mm(:,1)+offset,coords_mm(:,2)-offset,coords_mm(:,3)];
            merstruct.currentmer.posterior.trajectory{side} = [trajectory(:,1)+offset,trajectory(:,2)-offset,trajectory(:,3)];
            merstruct.currentmer.lateral.coords_mm{side} =  [coords_mm(:,1)+offset*2,coords_mm(:,2),coords_mm(:,3)]; %2mm medial (right is positive)
            merstruct.currentmer.lateral.trajectory{side} = [trajectory(:,1)+offset*2,trajectory(:,2),trajectory(:,3)]; %2mm medial (right is positive)
        elseif side==2 %left
            merstruct.currentmer.central.coords_mm{side} = [coords_mm(:,1)-offset,coords_mm(:,2),coords_mm(:,3)];
            merstruct.currentmer.central.trajectory{side} = [trajectory(:,1)-offset,trajectory(:,2),trajectory(:,3)];
            merstruct.currentmer.anterior.coords_mm{side} = [coords_mm(:,1)-offset,coords_mm(:,2)+offset,coords_mm(:,3)];
            merstruct.currentmer.anterior.trajectory{side} = [trajectory(:,1)-offset,trajectory(:,2)+offset,trajectory(:,3)];
            merstruct.currentmer.posterior.coords_mm{side} = [coords_mm(:,1)-offset,coords_mm(:,2)-offset,coords_mm(:,3)];
            merstruct.currentmer.posterior.trajectory{side} = [trajectory(:,1)-offset,trajectory(:,2)-offset,trajectory(:,3)];
            merstruct.currentmer.lateral.coords_mm{side} =  [coords_mm(:,1)-offset*2,coords_mm(:,2),coords_mm(:,3)]; %2mm medial (right is positive)
            merstruct.currentmer.lateral.trajectory{side} = [trajectory(:,1)-offset*2,trajectory(:,2),trajectory(:,3)]; %2mm medial (right is positive)
        end
end


function mer = ea_coords2mer(coords_mm,trajectory,resultfig,options)
        
offset = getfield(getappdata(resultfig,'merstruct'),'offset');
contact_length = getfield(getappdata(resultfig,'elspec'),'contact_length');

% x-axis --> negative = medial, positive = lateral
% y-axis --> negative = posterior, positive = anterior
% z-axis --> negative = inferior, positive = superior
for side=1:length(options.sides)
    dxyz = sqrt((diff(coords_mm{side}(1:2,1))^2)+(diff(coords_mm{side}(1:2,2))^2)+diff(coords_mm{side}(1:2,3))^2);
    slope = mean(diff(coords_mm{side}))/dxyz; %mean(diff(coords_mm{side}))/norm(mean(diff(coords_mm{side})))
    coords_mm{side} = coords_mm{side}-repmat(slope*contact_length/2,length(coords_mm{side}),1);
    trajectory{side} = trajectory{side}-repmat(slope*contact_length/2,length(trajectory{side}),1);
    if ea_getnativemni==1
        mer.central.coords_mm{side} = coords_mm{side};
        mer.central.trajectory{side} = trajectory{side};
        mer.anterior.coords_mm{side} = [coords_mm{side}(:,1),coords_mm{side}(:,2)+offset,coords_mm{side}(:,3)]; %2mm anterior
        mer.anterior.trajectory{side} = [trajectory{side}(:,1),trajectory{side}(:,2)+offset,trajectory{side}(:,3)]; %2mm anterior
        mer.posterior.coords_mm{side} = [coords_mm{side}(:,1),coords_mm{side}(:,2)-offset,coords_mm{side}(:,3)]; %2mm posterior
        mer.posterior.trajectory{side} = [trajectory{side}(:,1),trajectory{side}(:,2)-offset,trajectory{side}(:,3)]; %2mm posterior
        if side==1 %right
            mer.lateral.coords_mm{side} = [coords_mm{side}(:,1)+offset,coords_mm{side}(:,2),coords_mm{side}(:,3)]; %2mm lateral (right is positive)
            mer.lateral.trajectory{side} = [trajectory{side}(:,1)+offset,trajectory{side}(:,2),trajectory{side}(:,3)]; %2mm lateral (right is positive)
            mer.medial.coords_mm{side} = [coords_mm{side}(:,1)-offset,coords_mm{side}(:,2),coords_mm{side}(:,3)]; %2mm medial (right is positive)
            mer.medial.trajectory{side} = [trajectory{side}(:,1)-offset,trajectory{side}(:,2),trajectory{side}(:,3)]; %2mm medial (right is positive)
        elseif side==2 %left
            mer.lateral.coords_mm{side} = [coords_mm{side}(:,1)-offset,coords_mm{side}(:,2),coords_mm{side}(:,3)]; %2mm lateral (left is negative)
            mer.lateral.trajectory{side} = [trajectory{side}(:,1)-offset,trajectory{side}(:,2),trajectory{side}(:,3)]; %2mm lateral (left is negative)
            mer.medial.coords_mm{side} = [coords_mm{side}(:,1)+offset,coords_mm{side}(:,2),coords_mm{side}(:,3)]; %2mm medial (left is negative)
            mer.medial.trajectory{side} = [trajectory{side}(:,1)+offset,trajectory{side}(:,2),trajectory{side}(:,3)]; %2mm medial (left is negative)
        end
    end
    
    if ea_getnativemni==2
        mer.central.coords_mm{side} = coords_mm{side};
        mer.central.trajectory{side} = trajectory{side};
        mer.anterior.coords_mm{side} = [coords_mm{side}(:,1),coords_mm{side}(:,2)+offset,coords_mm{side}(:,3)]; %2mm anterior
        mer.anterior.trajectory{side} = [trajectory{side}(:,1),trajectory{side}(:,2)+offset,trajectory{side}(:,3)]; %2mm anterior
        mer.posterior.coords_mm{side} = [coords_mm{side}(:,1),coords_mm{side}(:,2)-offset,coords_mm{side}(:,3)]; %2mm posterior
        mer.posterior.trajectory{side} = [trajectory{side}(:,1),trajectory{side}(:,2)-offset,trajectory{side}(:,3)]; %2mm posterior
        if side==1 %right
            mer.lateral.coords_mm{side} = [coords_mm{side}(:,1)+offset,coords_mm{side}(:,2),coords_mm{side}(:,3)]; %2mm lateral (right is positive)
            mer.lateral.trajectory{side} = [trajectory{side}(:,1)+offset,trajectory{side}(:,2),trajectory{side}(:,3)]; %2mm lateral (right is positive)
            mer.medial.coords_mm{side} = [coords_mm{side}(:,1)-offset,coords_mm{side}(:,2),coords_mm{side}(:,3)]; %2mm medial (right is positive)
            mer.medial.trajectory{side} = [trajectory{side}(:,1)-offset,trajectory{side}(:,2),trajectory{side}(:,3)]; %2mm medial (right is positive)
        elseif side==2 %left
            mer.lateral.coords_mm{side} = [coords_mm{side}(:,1)-offset,coords_mm{side}(:,2),coords_mm{side}(:,3)]; %2mm lateral (left is negative)
            mer.lateral.trajectory{side} = [trajectory{side}(:,1)-offset,trajectory{side}(:,2),trajectory{side}(:,3)]; %2mm lateral (left is negative)
            mer.medial.coords_mm{side} = [coords_mm{side}(:,1)+offset,coords_mm{side}(:,2),coords_mm{side}(:,3)]; %2mm medial (left is negative)
            mer.medial.trajectory{side} = [trajectory{side}(:,1)+offset,trajectory{side}(:,2),trajectory{side}(:,3)]; %2mm medial (left is negative)
        end
    end
end
    


function refreshresultfig(handles,resultfig,options,side)
% this part makes changes of the figure active:
clearmertrajectories(handles,resultfig,options,side)
merstruct=getappdata(resultfig,'merstruct');
merhandles=getappdata(resultfig,'merhandles');
mermarkers=getappdata(resultfig,'mermarkers');

% Draw default trajectories
set(0,'CurrentFigure',resultfig); hold on;
for iSide = side
    if iSide==1 % right = 1 
        sidestr=lower('right');
    elseif iSide==2 % left = 2
        sidestr=lower('left');
    end
    for pos_ix = 1:length(merstruct.tract_labels)
        pos = merstruct.tract_labels{pos_ix};
        traj = merstruct.currentmer.(pos).trajectory{side};
        merhandles.(pos){side} = plot3(traj(:,1), traj(:,2), traj(:,3),...
            'color',merstruct.colormap(pos_ix,:), 'linew',5, 'tag',[pos '_' sidestr]);
    end
    set(handles.editimplanteddepth_left,'String','0');
end

setappdata(resultfig,'merhandles',merhandles)


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


function ea_keypress(handles,event)
% this is the main keypress function for the resultfigure. Add event
% listeners here.

if ismember('shift',event.Modifier)
end


function ea_update_trajectory_helper(hObject, handles, side, pos_str)
if side==1
    side_str = 'right';
else
    side_str = 'left';
end
resultfig=getappdata(handles.mercontrolfig,'resultfig');
merstruct=getappdata(resultfig,'merstruct');
dist = str2double(get(hObject,'String'))-str2double(get(handles.(['editimplanteddepth_' side_str]),'String'));
trajectory = ea_getmertrajectory(merstruct.currentmer.(pos_str).coords_mm{side},dist,merstruct.length,50);
ea_updatemertrajectory(handles,trajectory,0,['pos' pos_str '_' side_str]);


function ea_set_implanted_helper(handles, side)
if side == 1
    side_str = 'right';
else
    side_str = 'left';
end
resultfig=getappdata(handles.mercontrolfig,'resultfig');
merstruct=getappdata(resultfig,'merstruct');
options=getappdata(handles.mercontrolfig,'options');
tractstring=lower(get(handles.(['popupimplantedtract_' side_str]),'String'));
tracttag=tractstring{handles.(['popupimplantedtract_' side_str]).Value};
merstruct = ea_updatetracts(merstruct,side,tracttag,resultfig,handles);
for pos_str = merstruct.tract_labels
    ea_updatemertrajectory(handles,...
        merstruct.currentmer.(pos_str{:}).trajectory{side},...
        [pos_str{:} '_' side_str]);
end
setappdata(resultfig,'merstruct',merstruct);


% --- Executes on button press in pushbutton_clear_checks.
function pushbutton_clear_checks_Callback(hObject, eventdata, handles)
clearkeyboardcontrolgroup(handles);


function clearkeyboardcontrolgroup(handles)
resultfig = getappdata(handles.mercontrolfig, 'resultfig');
merstruct = getappdata(resultfig,'merstruct');
ui_tags = {handles.keyboardcontrolgroup.Children.Tag};
for side_str = {'left', 'right'}
    for tract_info = merstruct.tract_info
        set(handles.keyboardcontrolgroup.Children(strcmpi(ui_tags,...
                    ['keycheck' tract_info.label '_' side_str{:}])),...
            'Value', 0);
    end
end
setappdata(resultfig,'keymer',[])
