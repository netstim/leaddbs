function varargout = ea_select_colormap(varargin)
% EA_SELECT_COLORMAP MATLAB code for ea_select_colormap.fig
%      EA_SELECT_COLORMAP, by itself, creates a new EA_SELECT_COLORMAP or raises the existing
%      singleton*.
%
%      H = EA_SELECT_COLORMAP returns the handle to a new EA_SELECT_COLORMAP or the handle to
%      the existing singleton*.
%
%      EA_SELECT_COLORMAP('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in EA_SELECT_COLORMAP.M with the given input arguments.
%
%      EA_SELECT_COLORMAP('Property','Value',...) creates a new EA_SELECT_COLORMAP or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ea_select_colormap_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ea_select_colormap_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ea_select_colormap

% Last Modified by GUIDE v2.5 23-Mar-2020 23:05:53

% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ea_select_colormap_OpeningFcn, ...
                   'gui_OutputFcn',  @ea_select_colormap_OutputFcn, ...
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


% --- Executes just before ea_select_colormap is made visible.
function ea_select_colormap_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ea_select_colormap (see VARARGIN)

% Choose default command line output for ea_select_colormap
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

presets = {'parula'
           'jet'
           'hsv'
           'hot'
           'cool'
           'spring'
           'summer'
           'autumn'
           'winter'
           'gray'
           'bone'
           'copper'
           'pink'};
set(handles.presetcolormap, 'String', presets);

setappdata(handles.custom1, 'custom1startcolor', [1,1,1]);   % White
setappdata(handles.custom1, 'custom1endcolor',  [1,0,0]);    % Red
setappdata(handles.custom2, 'custom2startcolor',  [0,0,1]);  % Blue
setappdata(handles.custom2, 'custom2middlecolor',  [1,1,1]); % White
setappdata(handles.custom2, 'custom2endcolor',  [1,0,0]);    % Red

if ~isempty(varargin)
    if ~isempty(varargin{1})
        set(handles.colormapsize, 'String', num2str(varargin{1}));
    end

    if length(varargin) >= 2
        switch(varargin{2})
            case presets
                set(handles.colormapsetting, 'SelectedObject', handles.preset);
                set(handles.presetcolormap, 'Enable', 'on');
                set(handles.custom1startcolor, 'Enable', 'off');
                set(handles.custom1endcolor, 'Enable', 'off');
                set(handles.custom2startcolor, 'Enable', 'off');
                set(handles.custom2middlecolor, 'Enable', 'off');
                set(handles.custom2endcolor, 'Enable', 'off');

                set(handles.presetcolormap, 'Value', find(ismember(presets,varargin{2})));
                cmap = eval([varargin{2}, '(', get(handles.colormapsize, 'String') ,')']);
            case 'custom1'
                set(handles.colormapsetting, 'SelectedObject', handles.custom1);
                set(handles.presetcolormap, 'Enable', 'off');
                set(handles.custom1startcolor, 'Enable', 'on');
                set(handles.custom1endcolor, 'Enable', 'on');
                set(handles.custom2startcolor, 'Enable', 'off');
                set(handles.custom2middlecolor, 'Enable', 'off');
                set(handles.custom2endcolor, 'Enable', 'off');

                if length(varargin) == 3
                    cmapsize = length(varargin{3});
                    set(handles.colormapsize, 'String', num2str(cmapsize));
                    cmap = varargin{3};
                else
                    color1 = getappdata(handles.custom1, 'custom1startcolor');
                    color2 = getappdata(handles.custom1, 'custom1endcolor');

                    cmapsize = str2double(get(handles.colormapsize, 'String'));
                    cmap = ea_colorgradient(cmapsize, color1, color2);
                end

                setappdata(handles.selectcolormap, 'colormap', cmap);
            case 'custom2'
                set(handles.colormapsetting, 'SelectedObject', handles.custom2);
                set(handles.presetcolormap, 'Enable', 'off');
                set(handles.custom1startcolor, 'Enable', 'off');
                set(handles.custom1endcolor, 'Enable', 'off');
                set(handles.custom2startcolor, 'Enable', 'on');
                set(handles.custom2middlecolor, 'Enable', 'on');
                set(handles.custom2endcolor, 'Enable', 'on');

                if length(varargin) == 3
                    cmapsize = length(varargin{3});
                    set(handles.colormapsize, 'String', num2str(cmapsize));
                    cmap = varargin{3};
                else
                    color1 = getappdata(handles.custom2, 'custom2startcolor');
                    color2 = getappdata(handles.custom2, 'custom2middlecolor');
                    color3 = getappdata(handles.custom2, 'custom2endcolor');
                    cmapsize = str2double(get(handles.colormapsize, 'String'));
                    cmap = ea_colorgradient(cmapsize, color1, color2, color3);
                end
                setappdata(handles.selectcolormap, 'colormap', cmap);
        end
    end
end

if length(varargin) < 2
    cmap = presets{get(handles.presetcolormap, 'Value')};
    cmap = eval([cmap, '(', get(handles.colormapsize, 'String') ,')']);
end

setappdata(handles.selectcolormap, 'colormap', cmap);
updatePreview(cmap, handles.previewcolorbar);

% UIWAIT makes ea_select_colormap wait for user response (see UIRESUME)
uiwait(handles.selectcolormap);


% --- Executes when user attempts to close selectcolormap.
function selectcolormap_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to selectcolormap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Choose default command line output for ea_select_colormap
handles.output = getappdata(handles.selectcolormap, 'colormap');

% Update handles structure
guidata(hObject, handles);

uiresume(handles.selectcolormap);


% --- Outputs from this function are returned to the command line.
function varargout = ea_select_colormap_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure

varargout{1} = handles.output;
delete(hObject);


% --- Executes during object creation, after setting all properties.
function colormapsize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to colormapsize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
set(hObject, 'String', num2str(length(gray)));


function colormapsize_Callback(hObject, eventdata, handles)
% hObject    handle to colormapsize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of colormapsize as text
%        str2double(get(hObject,'String')) returns contents of colormapsize as a double
switch get(get(handles.colormapsetting, 'SelectedObject'), 'Tag');
    case 'preset'
        presets = get(handles.presetcolormap, 'String');
        cmap = presets{get(handles.presetcolormap, 'Value')};
        cmap = eval([cmap, '(', get(handles.colormapsize, 'String') ,')']);
    case 'custom1'
        color1 = getappdata(handles.custom1, 'custom1startcolor');
        color2 = getappdata(handles.custom1, 'custom1endcolor');

        cmapsize = str2double(get(handles.colormapsize, 'String'));
        cmap = ea_colorgradient(cmapsize, color1, color2);
    case 'custom2'
        color1 = getappdata(handles.custom2, 'custom2startcolor');
        color2 = getappdata(handles.custom2, 'custom2middlecolor');
        color3 = getappdata(handles.custom2, 'custom2endcolor');

        cmapsize = str2double(get(handles.colormapsize, 'String'));
        cmap = ea_colorgradient(cmapsize, color1, color2, color3);
end

setappdata(handles.selectcolormap, 'colormap', cmap);
updatePreview(cmap, handles.previewcolorbar);


% --- Executes on button press in preset.
function preset_Callback(hObject, eventdata, handles)
% hObject    handle to preset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of preset
set(handles.presetcolormap, 'Enable', 'on');
set(handles.custom1startcolor, 'Enable', 'off');
set(handles.custom1endcolor, 'Enable', 'off');
set(handles.custom2startcolor, 'Enable', 'off');
set(handles.custom2middlecolor, 'Enable', 'off');
set(handles.custom2endcolor, 'Enable', 'off');

presets = get(handles.presetcolormap, 'String');
cmap = presets{get(handles.presetcolormap, 'Value')};
cmap = eval([cmap, '(', get(handles.colormapsize, 'String') ,')']);
setappdata(handles.selectcolormap, 'colormap', cmap);
updatePreview(cmap, handles.previewcolorbar);


% --- Executes during object creation, after setting all properties.
function presetcolormap_CreateFcn(hObject, eventdata, handles)
% hObject    handle to presetcolormap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in presetcolormap.
function presetcolormap_Callback(hObject, eventdata, handles)
% hObject    handle to presetcolormap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns presetcolormap contents as cell array
%        contents{get(hObject,'Value')} returns selected item from presetcolormap
presets = get(handles.presetcolormap, 'String');
cmap = presets{get(handles.presetcolormap, 'Value')};
cmap = eval([cmap, '(', get(handles.colormapsize, 'String') ,')']);
setappdata(handles.selectcolormap, 'colormap', cmap);
updatePreview(cmap, handles.previewcolorbar);


% --- Executes on button press in custom1.
function custom1_Callback(hObject, eventdata, handles)
% hObject    handle to custom1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of custom1
set(handles.presetcolormap, 'Enable', 'off');
set(handles.custom1startcolor, 'Enable', 'on');
set(handles.custom1endcolor, 'Enable', 'on');
set(handles.custom2startcolor, 'Enable', 'off');
set(handles.custom2middlecolor, 'Enable', 'off');
set(handles.custom2endcolor, 'Enable', 'off');

color1 = getappdata(handles.custom1, 'custom1startcolor');
color2 = getappdata(handles.custom1, 'custom1endcolor');

cmapsize = str2double(get(handles.colormapsize, 'String'));
cmap = ea_colorgradient(cmapsize, color1, color2);
setappdata(handles.selectcolormap, 'colormap', cmap);
updatePreview(cmap, handles.previewcolorbar);
guidata(hObject,handles)


% --- Executes on button press in custom1startcolor.
function custom1startcolor_Callback(hObject, eventdata, handles)
% hObject    handle to custom1startcolor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
setappdata(handles.custom1, 'custom1startcolor',  ea_uisetcolor([1,1,1]));
cmapsize = str2double(get(handles.colormapsize, 'String'));
color1 = getappdata(handles.custom1, 'custom1startcolor');
color2 = getappdata(handles.custom1, 'custom1endcolor');
cmap = ea_colorgradient(cmapsize, color1, color2);
setappdata(handles.selectcolormap, 'colormap', cmap);
updatePreview(cmap, handles.previewcolorbar);


% --- Executes on button press in custom1endcolor.
function custom1endcolor_Callback(hObject, eventdata, handles)
% hObject    handle to custom1endcolor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
setappdata(handles.custom1, 'custom1endcolor',  ea_uisetcolor([1,0,0]));

cmapsize = str2double(get(handles.colormapsize, 'String'));
color1 = getappdata(handles.custom1, 'custom1startcolor');
color2 = getappdata(handles.custom1, 'custom1endcolor');
cmap = ea_colorgradient(cmapsize, color1, color2);
setappdata(handles.selectcolormap, 'colormap', cmap);
updatePreview(cmap, handles.previewcolorbar);


% --- Executes on button press in custom1.
function custom2_Callback(hObject, eventdata, handles)
% hObject    handle to custom2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of custom1
set(handles.presetcolormap, 'Enable', 'off');
set(handles.custom1startcolor, 'Enable', 'off');
set(handles.custom1endcolor, 'Enable', 'off');
set(handles.custom2startcolor, 'Enable', 'on');
set(handles.custom2middlecolor, 'Enable', 'on');
set(handles.custom2endcolor, 'Enable', 'on');

color1 = getappdata(handles.custom2, 'custom2startcolor');
color2 = getappdata(handles.custom2, 'custom2middlecolor');
color3 = getappdata(handles.custom2, 'custom2endcolor');

cmapsize = str2double(get(handles.colormapsize, 'String'));
cmap = ea_colorgradient(cmapsize, color1, color2, color3);
setappdata(handles.selectcolormap, 'colormap', cmap);
updatePreview(cmap, handles.previewcolorbar);


% --- Executes on button press in custom2startcolor.
function custom2startcolor_Callback(hObject, eventdata, handles)
% hObject    handle to custom2startcolor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
setappdata(handles.custom2, 'custom2startcolor',  ea_uisetcolor([0,0,1]));

cmapsize = str2double(get(handles.colormapsize, 'String'));
color1 = getappdata(handles.custom2, 'custom2startcolor');
color2 = getappdata(handles.custom2, 'custom2middlecolor');
color3 = getappdata(handles.custom2, 'custom2endcolor');
cmap = ea_colorgradient(cmapsize, color1, color2, color3);
setappdata(handles.selectcolormap, 'colormap', cmap);
updatePreview(cmap, handles.previewcolorbar);


% --- Executes on button press in custom2middlecolor.
function custom2middlecolor_Callback(hObject, eventdata, handles)
% hObject    handle to custom2middlecolor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
setappdata(handles.custom2, 'custom2middlecolor',  ea_uisetcolor([1,1,1]));

cmapsize = str2double(get(handles.colormapsize, 'String'));
color1 = getappdata(handles.custom2, 'custom2startcolor');
color2 = getappdata(handles.custom2, 'custom2middlecolor');
color3 = getappdata(handles.custom2, 'custom2endcolor');
cmap = ea_colorgradient(cmapsize, color1, color2, color3);
setappdata(handles.selectcolormap, 'colormap', cmap);
updatePreview(cmap, handles.previewcolorbar);


% --- Executes on button press in custom2endcolor.
function custom2endcolor_Callback(hObject, eventdata, handles)
% hObject    handle to custom2endcolor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
setappdata(handles.custom2, 'custom2endcolor',  ea_uisetcolor([1,0,0]));

cmapsize = str2double(get(handles.colormapsize, 'String'));
color1 = getappdata(handles.custom2, 'custom2startcolor');
color2 = getappdata(handles.custom2, 'custom2middlecolor');
color3 = getappdata(handles.custom2, 'custom2endcolor');
cmap = ea_colorgradient(cmapsize, color1, color2, color3);
setappdata(handles.selectcolormap, 'colormap', cmap);
updatePreview(cmap, handles.previewcolorbar);


% --- Executes on button press in confirm.
function confirm_Callback(hObject, eventdata, handles)
% hObject    handle to confirm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close(handles.selectcolormap);


function updatePreview(cmap, axes)
map = colormap(axes, cmap);
width = ceil(length(cmap)/16);
image(repmat(cat(3, map(:,1)', map(:,2)', map(:,3)'), width, 1));

% Remove yticks
set(gca, 'ytick', []);
set(gca, 'xtick', []);
