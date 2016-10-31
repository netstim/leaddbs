function varargout = mrLoadRetGUI(varargin)
% fig = mrLoadRetGUI('viewNum',viewNum)
% $Id: mrLoadRetGUI.m 2962 2014-09-11 22:04:32Z justin $
%
% Creates a new mrLoadRet GUI.
% This function was created along with mrLoadRetGui.fig using GUIDE.
%
% djh, 6/2004
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Last Modified by GUIDE v2.5 01-Aug-2012 14:25:01
% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @mrLoadRetGUI_OpeningFcn, ...
    'gui_OutputFcn',  @mrLoadRetGUI_OutputFcn, ...
    'gui_LayoutFcn',  [], ...
    'gui_Callback',   []);
if nargin & isstr(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

% --- Executes just before mrLoadRetGUI is made visible.
function mrLoadRetGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

% Parse varargin
% mrLoadRetGUI must be called as follows:
%    mrLoadRetGUI('viewNum',viewNum)
% The view structure must already exist in MRL.views{viewNum}
for index = 1:2:length(varargin)
    field = varargin{index};
    val = varargin{index+1};
    switch field
        case 'viewNum'
            viewNum = val;
            handles.viewNum = viewNum;
        otherwise
            warn('mrLoadRetGUI: invalid initialization argument')
    end
end

% Initialize coords and slice orientation
handles.coords = [1,1,1];
handles.sliceOrientation = 3;
guidata(hObject,handles);

% Initialize group popup
set(handles.groupPopup,'String',viewGet([],'groupNames'));

% Initialize rotate slider
set(handles.rotateSlider,'sliderStep',[1/360 10/360]);

% Enable/disable various widgets depending on viewType
% removed the switch here, since we always have one viewType

% Initialize the slice orientation radio buttons
set(handles.sagittalRadioButton,'Value',1);
set(handles.coronalRadioButton,'Value',0);
set(handles.axialRadioButton,'Value',0);
set(handles.corticalDepth,'Visible','off');
set(handles.corticalDepthSlider,'Visible','off');
set(handles.corticalDepthText,'Visible','off');
if isfield(handles,'corticalMaxDepthSlider')
  set(handles.corticalMaxDepthSlider,'Visible','off');
  set(handles.corticalMaxDepthText,'Visible','off');
end

% Choose default command line output for mrLoadRetGUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% --------------------------------------------------------------------
% Output function
function varargout = mrLoadRetGUI_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);

% Get default command line output from handles structure
varargout{1} = handles.output;

% --------------------------------------------------------------------
% Button Down functions

function axis_ButtonDownFcn(hObject, eventdata, handles)

function figure_ButtonDownFcn(hObject, eventdata, handles)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    windowButtonMotion function    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function figure_WindowButtonMotionFcn(hObject, eventdata, handles)

gui = guidata(hObject);
v = viewGet([],'view',gui.viewNum);
if ~isempty(viewGet(v,'curBase')) && isequal(1,viewGet(v,'baseMultiAxis'))
  % then get coords
  coords = mlrGetMouseCoords(v);
  if (coords.inAxis ~= -1)
    set(hObject,'Pointer','crosshair');
  else
    set(hObject,'Pointer','arrow');
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    windowButtonDown function    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function figure_WindowButtonDownFcn(hObject, eventdata, handles)

% get the view
v = viewGet(handles.viewNum,'view');

% see if we are displaying a base and are displaying all planes
baseType = viewGet(v,'baseType');
if ~isempty(baseType) && isequal(viewGet(v,'baseType'),0) && isequal(viewGet(v,'baseMultiAxis'),1)
  % then get coords
  coords = mlrGetMouseCoords(v);
  if ~isempty(coords.base)
    % set the change in coords
    curCoords = [coords.base.x coords.base.y coords.base.z];
    v = viewSet(v,'curCoords',curCoords);
    v = viewSet(v,'curSlice',curCoords(viewGet(v,'baseSliceIndex')));
    % and display
    refreshMLRDisplay(viewGet(v,'viewNum'));
  end
end

%%%%%%%%%%%%%%%%
%    Resize    %
%%%%%%%%%%%%%%%%
function figure_ResizeFcn(hObject, eventdata, handles)
% Change the axis size to fill the figure, leaving room at the bottom for
% the widgets.
if exist('handles','var') & ~isempty(handles)
    figureSize = get(handles.figure,'Position');
    axisSize = get(handles.axis,'Position');
    axisSize(3) = figureSize(3);
    axisSize(4) = figureSize(4) - axisSize(2);
    set(handles.axis,'Position',axisSize);
    refreshMLRDisplay(viewNum);
end


%%%%%%%%%%%%%%%%%%%%%%%%%
%    Create function    %
%%%%%%%%%%%%%%%%%%%%%%%%%
function figure_CreateFcn(hObject, eventdata, handles)

% hook the above button motion function in
set(hObject,'WindowButtonMotionFcn',@figure_WindowButtonMotionFcn);

% --------------------------------------------------------------------
% Delete  and Close functions

function axis_DeleteFcn(hObject, eventdata, handles)

function figure_DeleteFcn(hObject, eventdata, handles)
if ~ieNotDefined('handles') & isfield(handles,'viewNum')
    % Delete the view
    viewNum = handles.viewNum;
    deleteView(viewNum);
end

function figure_CloseRequestFcn(hObject, eventdata, handles)
% Hint: delete(hObject) closes the figure
delete(hObject);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Radio buttons

% --- Sagittal
function sagittalRadioButton_Callback(hObject, eventdata, handles)
viewNum = handles.viewNum;
viewSet(viewNum,'sliceOrientation','sagittal');
if ~viewGet(viewNum,'baseMultiAxis') %only in single slice view
  refreshMLRDisplay(viewNum);
end

% --- Coronal
function coronalRadioButton_Callback(hObject, eventdata, handles)
viewNum = handles.viewNum;
viewSet(viewNum,'sliceOrientation','coronal');
if ~viewGet(viewNum,'baseMultiAxis') %only in single slice view
  refreshMLRDisplay(viewNum);
end

% --- Axial
function axialRadioButton_Callback(hObject, eventdata, handles)
viewNum = handles.viewNum;
viewSet(viewNum,'sliceOrientation','axial');
if ~viewGet(viewNum,'baseMultiAxis') %only in single slice view
  refreshMLRDisplay(viewNum);
end

% --- Left
function leftRadioButton_Callback(hObject, eventdata, handles)
% *** Not tested ***
viewNum = handles.viewNum;
mlrGuiSet(viewNum,'slice',1);
refreshMLRDisplay(viewNum);

% --- Right
function rightRadioButton_Callback(hObject, eventdata, handles)
% *** Not tested ***
viewNum = handles.viewNum;
mlrGuiSet(viewNum,'slice',2);
refreshMLRDisplay(viewNum);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Popups

% --- Group
function groupPopup_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function groupPopup_Callback(hObject, eventdata, handles)
mrGlobals
viewNum = handles.viewNum;
% Update the group number
viewSet(viewNum,'curGroup',get(hObject,'Value'));
% Delete the overlays
%MLR.views{viewNum}.analyses = [];
%MLR.views{viewNum}.curAnalysis = [];
% Update nScans
mlrGuiSet(viewNum,'nScans',viewGet(viewNum,'nScans'));
% Refresh the display
refreshMLRDisplay(viewNum);

% --- ROI
function roiPopup_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function roiPopup_Callback(hObject, eventdata, handles)
viewNum = handles.viewNum;
viewSet(viewNum,'curROI',get(hObject,'Value'));
refreshMLRDisplay(viewNum);

% --- Analysis
function analysisPopup_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function analysisPopup_Callback(hObject, eventdata, handles)
viewNum = handles.viewNum;
viewSet(viewNum,'curAnalysis',get(hObject,'Value'));
refreshMLRDisplay(viewNum);

% --- Overlay
function overlayPopup_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function overlayPopup_Callback(hObject, eventdata, handles)
viewNum = handles.viewNum;
viewSet(viewNum,'curOverlay',get(hObject,'Value'));
refreshMLRDisplay(viewNum);

% --- Base image
function basePopup_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function basePopup_Callback(hObject, eventdata, handles)
viewNum = handles.viewNum;
viewSet(viewNum,'curBase',get(hObject,'Value'));
refreshMLRDisplay(viewNum);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sliders and associated text fields

% --- Scan
function scanSlider_CreateFcn(hObject, eventdata, handles)
usewhitebg = 1;
if usewhitebg
    set(hObject,'BackgroundColor',[.9 .9 .9]);
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function scanText_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function scanSlider_Callback(hObject, eventdata, handles)

% get view/viewNum
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};

value = round(get(hObject,'Value'));
mlrGuiSet(viewNum,'scanText',value);
view = viewSet(view,'curScan',value);
refreshMLRDisplay(viewNum);

function scanText_Callback(hObject, eventdata, handles)

% get view/viewNum
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};

value = round(str2num(get(hObject,'String')));
if length(value)~=1 %if the user just erased the value, get it from the view and do nothing
  set(hObject,'String',num2str(viewGet(view,'curScan')));
else %otherwise, set the new value in the view and the GUI
  %set the current scan in the view
  view = viewSet(view,'curScan',value);
  %set the current scan on slider
  mlrGuiSet(viewNum,'scan',viewGet(view,'curScan'));
  refreshMLRDisplay(viewNum);
end

% --- Slice
function sliceSlider_CreateFcn(hObject, eventdata, handles)
usewhitebg = 1;
if usewhitebg
    set(hObject,'BackgroundColor',[.9 .9 .9]);
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function sliceText_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function sliceSlider_Callback(hObject, eventdata, handles)

% get view/viewNum
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};

value = round(get(hObject,'Value'));
mlrGuiSet(viewNum,'sliceText',value); %set slice edit box value
viewSet(view,'curSlice',value); % set the current slice
refreshMLRDisplay(viewNum);

function sliceText_Callback(hObject, eventdata, handles)

% get view/viewNum
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};

value = round(str2num(get(hObject,'String')));
if length(value)~=1 %if the user just erased the value, get it from the slider and do nothing
  set(hObject,'String',num2str(get(handles.sliceSlider,'value')));
else %otherwise, set the new value in the GUI and then in the view
  mlrGuiSet(viewNum,'slice',value); %doing this first ensures that the value is not outside the values permitted by the slider
  viewSet(view,'curSlice',get(handles.sliceSlider,'Value')); %use the value from the slider since it might have been changed by mlrGuiSet
  refreshMLRDisplay(viewNum);
end

% --- Executes on slider movement.
function corticalDepthSlider_Callback(hObject, eventdata, handles)

viewNum = handles.viewNum;
newMinValue = get(hObject,'Value');
if isfield(handles,'linkMinMaxDepthCheck') && get(handles.linkMinMaxDepthCheck,'value')
  maxDepth = viewGet(viewNum,'corticalMaxDepth');
  minDepth = str2num(get(handles.corticalDepthText,'string'));
  newMaxValue = maxDepth+newMinValue-minDepth;
  if newMaxValue>=0 && newMaxValue<=1 %if both values are in [0 1]
    viewSet(viewNum,'corticalMinDepth',newMinValue);
    viewSet(viewNum,'corticalMaxDepth',newMaxValue);
    drawnow;
    refreshMLRDisplay(viewNum);
  else
    set(hObject,'Value',minDepth);
  end
elseif isfield(handles,'corticalMaxDepthSlider')
  viewSet(viewNum,'corticalMinDepth',newMinValue);
  refreshMLRDisplay(viewNum);
else
  viewSet(viewNum,'corticalDepth',newMinValue);
  refreshMLRDisplay(viewNum);
end

% --- Executes during object creation, after setting all properties.
function corticalDepthSlider_CreateFcn(hObject, eventdata, handles)

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
  set(hObject,'BackgroundColor',[.9 .9 .9]);
end
corticalDepthBins=mrGetPref('corticaldepthbins');
corticalDepth=round((corticalDepthBins-1)/2)/(corticalDepthBins-1);
set(hObject,'value',corticalDepth);
set(hObject,'SliderStep',min(1/(corticalDepthBins-1)*[1 3],1));

function corticalDepthText_Callback(hObject, eventdata, handles)

% Hints: get(hObject,'String') returns contents of corticalDepthText as text
%        str2double(get(hObject,'String')) returns contents of corticalDepthText as a double

viewNum = handles.viewNum;
newMinValue = str2num(get(hObject,'String'));
if length(newMinValue) ~= 1%if the user just erased the value, get it from the slider and do nothing
  set(hObject,'String',num2str(get(handles.corticalDepthSlider,'value')));
else %otherwise, set the new value in the GUI (there is no field for cortical depth in the view, although there probably should)
  if isfield(handles,'linkMinMaxDepthCheck') && get(handles.linkMinMaxDepthCheck,'value')
    maxDepth = viewGet(viewNum,'corticalMaxDepth');
    minDepth = get(handles.corticalDepthSlider,'value');
    newMaxValue = maxDepth+newMinValue-minDepth;
    if newMaxValue>=0 && newMaxValue<=1 %if both values are in [0 1]
      viewSet(viewNum,'corticalMinDepth',newMinValue);
      viewSet(viewNum,'corticalMaxDepth',newMaxValue);
      refreshMLRDisplay(viewNum);
    else
      set(hObject,'string',num2str(minDepth));
    end
  elseif isfield(handles,'corticalMaxDepthSlider')
    viewSet(viewNum,'corticalMinDepth',newMinValue);
    refreshMLRDisplay(viewNum);
  else
    viewSet(viewNum,'corticalDepth',newMinValue);
    refreshMLRDisplay(viewNum);
  end
end

% --- Executes during object creation, after setting all properties.
function corticalDepthText_CreateFcn(hObject, eventdata, handles)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
corticalDepthBins=mrGetPref('corticaldepthbins');
corticalDepth=round((corticalDepthBins-1)/2)/(corticalDepthBins-1);
set(hObject,'string',num2str(corticalDepth));



% --- baseGamma
function baseGammaSlider_CreateFcn(hObject, eventdata, handles)
usewhitebg = 1;
if usewhitebg
    set(hObject,'BackgroundColor',[.9 .9 .9]);
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function baseGammaText_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function baseGammaSlider_Callback(hObject, eventdata, handles)
viewNum = handles.viewNum;
value = get(hObject,'Value');
viewSet(handles.viewNum,'baseGamma',value);
refreshMLRDisplay(viewNum);

function baseGammaText_Callback(hObject, eventdata, handles)
viewNum = handles.viewNum;
value = str2num(get(hObject,'String'));
if length(value)~=1 %if the user just erased the value, get it from the slider and do nothing
  set(hObject,'String',num2str(get(handles.baseGammaSlider,'value')));
else %otherwise, set the new value in the view and the GUI
  viewSet(viewNum,'baseGamma',value);
  refreshMLRDisplay(viewNum);
end

% --- baseTilt
function baseTiltSlider_CreateFcn(hObject, eventdata, handles)
usewhitebg = 1;
if usewhitebg
    set(hObject,'BackgroundColor',[.9 .9 .9]);
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function baseTiltText_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function baseTiltSlider_Callback(hObject, eventdata, handles)
viewNum = handles.viewNum;
value = get(hObject,'Value');
viewSet(viewNum,'tilt',value);
v = viewGet([],'view',viewNum);

if (viewGet(v,'baseType') == 2)
  setMLRViewAngle(v);
else
  refreshMLRDisplay(viewNum);
end

function baseTiltText_Callback(hObject, eventdata, handles)
viewNum = handles.viewNum;
value = str2num(get(hObject,'String'));
if length(value)~=1 %if the user just erased the value, get it from the slider and do nothing
  set(hObject,'String',num2str(get(handles.baseTiltSlider,'value')));
else %otherwise, set the new value in the view and the GUI
  if value < 0,value = 0;end
  if value > 360,value = 360;end
  viewSet(viewNum,'tilt',value);
  v = viewGet([],'view',viewNum);
  fig = viewGet(v,'figNum');
  gui = guidata(fig);

  if (viewGet(v,'baseType') == 2)
    setMLRViewAngle(v);
  else
    refreshMLRDisplay(viewNum);
  end
end

% --- overlayMax
function overlayMaxSlider_CreateFcn(hObject, eventdata, handles)
usewhitebg = 1;
if usewhitebg
    set(hObject,'BackgroundColor',[.9 .9 .9]);
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function overlayMaxText_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function overlayMaxSlider_Callback(hObject, eventdata, handles)
viewNum = handles.viewNum;
value = get(hObject,'Value');
viewSet(viewNum,'overlayMax',value,viewGet(viewNum,'curClippingOverlay'));
refreshMLRDisplay(viewNum);

function overlayMaxText_Callback(hObject, eventdata, handles)
viewNum = handles.viewNum;
value = str2num(get(hObject,'String'));
if length(value) ~= 1 %if the user just erased the value, get it from the slider and do nothing
  set(hObject,'String',num2str(get(handles.overlayMaxSlider,'value')));
else %otherwise, set the new value in the view and the GUI
  viewSet(viewNum,'overlayMax',value,viewGet(viewNum,'curClippingOverlay'));
  refreshMLRDisplay(viewNum);
end

% --- overlayMin
function overlayMinSlider_CreateFcn(hObject, eventdata, handles)
usewhitebg = 1;
if usewhitebg
    set(hObject,'BackgroundColor',[.9 .9 .9]);
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function overlayMinText_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function overlayMinSlider_Callback(hObject, eventdata, handles)
viewNum = handles.viewNum;
value = get(hObject,'Value');
viewSet(viewNum,'overlayMin',value,viewGet(viewNum,'curClippingOverlay'));
refreshMLRDisplay(viewNum);

function overlayMinText_Callback(hObject, eventdata, handles)
viewNum = handles.viewNum;
value = str2num(get(hObject,'String'));
if length(value) ~= 1%if the user just erased the value, get it from the slider and do nothing
  set(hObject,'String',num2str(get(handles.overlayMinSlider,'value')));
else %otherwise, set the new value in the view and the GUI
  viewSet(viewNum,'overlayMin',value,viewGet(viewNum,'curClippingOverlay'));
  refreshMLRDisplay(viewNum);
end

% --- alpha
function alphaSlider_CreateFcn(hObject, eventdata, handles)
usewhitebg = 1;
if usewhitebg
    set(hObject,'BackgroundColor',[.9 .9 .9]);
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function alphaText_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function alphaSlider_Callback(hObject, eventdata, handles)
viewNum = handles.viewNum;
value = get(hObject,'Value');
viewSet(viewNum,'alpha',value);
refreshMLRDisplay(viewNum);

function alphaText_Callback(hObject, eventdata, handles)
viewNum = handles.viewNum;
value = str2num(get(hObject,'String'));
if length(value)~=1 %if the user just erased the value, get it from the slider and do nothing
  set(hObject,'String',num2str(get(handles.alphaSlider,'value')));
else %otherwise, set the new value in the view and the GUI
  viewSet(viewNum,'alpha',value);
  refreshMLRDisplay(viewNum);
end

% --- rotate
function rotateSlider_CreateFcn(hObject, eventdata, handles)
usewhitebg = 1;
if usewhitebg
    set(hObject,'BackgroundColor',[.9 .9 .9]);
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function rotateText_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function rotateSlider_Callback(hObject, eventdata, handles)

viewNum = handles.viewNum;
value = get(hObject,'Value');
v = viewGet([],'view',viewNum);
v = viewSet(v,'rotate',value);

if (viewGet(v,'baseType') == 2)
  setMLRViewAngle(v);
else
  refreshMLRDisplay(viewNum);
end

function rotateText_Callback(hObject, eventdata, handles)
viewNum = handles.viewNum;
value = str2num(get(hObject,'String'));
if length(value)~=1 %if the user just erased the value, get it from the slider and do nothing
  set(hObject,'String',num2str(get(handles.rotateSlider,'value')));
else %otherwise, set the new value in the view and the GUI
  v = viewGet([],'view',viewNum);
  v = viewSet(v,'rotate',value);
  
  if (viewGet(v,'baseType') == 2)
    setMLRViewAngle(v);
  else
    refreshMLRDisplay(viewNum);
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fileMenu_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function fileBaseMenu_Callback(hObject, eventdata, handles)

function loadAnatomyMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = loadAnat(MLR.views{viewNum});
refreshMLRDisplay(viewNum);

% --------------------------------------------------------------------
function loadFromVolumeMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
% get the volume directory from prefs
volumeDirectory = mrGetPref('volumeDirectory');
saveVolumeDirectory = 0;
if isempty(volumeDirectory)
    saveVolumeDirectory = 1
end
% load the anatomy
[view volumeDirectory] = loadAnat(view,'',volumeDirectory);
refreshMLRDisplay(viewNum);
% if volumeDirectory prefMenu was empty before, than save it now
if saveVolumeDirectory
    mrSetPref('volumeDirectory',volumeDirectory);
end

% --------------------------------------------------------------------
function useCurrentScanMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
% get the view
viewNum = handles.viewNum;
v = MLR.views{viewNum};

% get the tseries path/filename
tSeriesPathStr = viewGet(v,'tSeriesPathStr',viewGet(v,'curScan'));
% load that as an anatome
if isfile(tSeriesPathStr)
    v = loadAnat(MLR.views{viewNum},getLastDir(tSeriesPathStr),fileparts(tSeriesPathStr));
    refreshMLRDisplay(viewNum);
else
    mrWarnDlg(sprintf('(mrLoadRetGUI) Could not find tSeries %s',tSeriesPathStr));
end

% --------------------------------------------------------------------
function importFlatMenuItem_Callback(hObject, eventdata, handles)

mrGlobals;
viewNum = handles.viewNum;
v = MLR.views{viewNum};
% v = loadFlat(v);
base = importFlatOFF;
if ~isempty(base)
  viewSet(v, 'newbase', base);
  refreshMLRDisplay(viewNum);
end

% --------------------------------------------------------------------
function importSurfaceMenuItem_Callback(hObject, eventdata, handles)

mrGlobals;
viewNum = handles.viewNum;
v = MLR.views{viewNum};
base = importSurfaceOFF;
if ~isempty(base)
  viewSet(v, 'newbase', base);
  refreshMLRDisplay(viewNum);
end

% --------------------------------------------------------------------
function saveAnatomyMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
n = viewGet(viewNum,'currentBase');
saveAnat(MLR.views{viewNum},n,1);

% --------------------------------------------------------------------
function SaveAsMenuItem_Callback(hObject, eventdata, handles)

mrGlobals;
viewNum = handles.viewNum;
n = viewGet(viewNum,'currentBase');
saveAnat(MLR.views{viewNum},n,1,1);

% --------------------------------------------------------------------
function fileAnalysisMenu_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function loadAnalysisMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = loadAnalysis(MLR.views{viewNum});
refreshMLRDisplay(viewNum);

% --------------------------------------------------------------------
function saveAnalysisMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
n = viewGet(viewNum,'currentAnalysis');
if ~isempty(n)
  saveAnalysis(MLR.views{viewNum},n);
end

% --------------------------------------------------------------------
function saveAllAnalysisMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
numberAnalyses = viewGet(viewNum,'numberofAnalyses');
for n = 1:numberAnalyses
    saveAnalysis(MLR.views{viewNum},n);
end

% --------------------------------------------------------------------
function fileOverlayMenu_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function loadOverlayMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = loadOverlay(MLR.views{viewNum});
refreshMLRDisplay(viewNum);

% --------------------------------------------------------------------
function saveOverlayMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
n = viewGet(viewNum,'currentOverlay');
m = viewGet(viewNum,'currentAnalysis');
saveOverlay(MLR.views{viewNum},n,m);

% --------------------------------------------------------------------
function saveAllOverlayMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
numberOverlays = viewGet(viewNum,'numberofOverlays');
m = viewGet(viewNum,'currentAnalysis');
for n = 1:numberOverlays
    saveOverlay(MLR.views{viewNum},n,m);
end

% --------------------------------------------------------------------
function fileRoiMenu_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function loadROIMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = loadROI(MLR.views{viewNum});
refreshMLRDisplay(viewNum);

% --------------------------------------------------------------------
function loadFromVolumeDirectoryROIMenuItem_Callback(hObject, eventdata, handles)

mrGlobals;
viewNum = handles.viewNum;
v = MLR.views{viewNum};
% get the volume directory from prefs
volumeDirectory = mrGetPref('volumeDirectory');
% load the rois
v = loadROI(v,[],[],volumeDirectory);
% and refresh
refreshMLRDisplay(viewNum);


% --------------------------------------------------------------------
function importROIMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
view = importROI(view);

% --------------------------------------------------------------------
function saveROIMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
for n = viewGet(viewNum,'currentroi');
  saveROI(MLR.views{viewNum},n,1);
end

% --------------------------------------------------------------------
function saveManyROIsMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
roiList = selectInList(MLR.views{viewNum},'rois','Select ROI(s) to save');
for iRoi = roiList
    saveROI(MLR.views{viewNum},iRoi,1);
end

% --------------------------------------------------------------------
function saveAllROIMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
numberROIs = viewGet(viewNum,'numberofrois');
for n = 1:numberROIs
    saveROI(MLR.views{viewNum},n,1);
end

% --------------------------------------------------------------------
function importMenu_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function importGroupMenuItem_Callback(hObject, eventdata, handles)

viewNum = handles.viewNum;
% turn off interrogator since importGroupScans has
% to temporarily switch the MLR views to get info from the import group
mrInterrogator('end',viewNum);
importGroupScans;

% --------------------------------------------------------------------
function importTSeriesMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
importTSeries(MLR.views{viewNum});

% --------------------------------------------------------------------
function importOverlayMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
importOverlay(MLR.views{viewNum});  %importOverlay checks the compatibilty of the imported data with the current scan

% --------------------------------------------------------------------
function exportMenu_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function exportAnatomyMenuItem_Callback(hObject, eventdata, handles)
mrWarnDlg('exportAnatomy not yet implemented');

% --------------------------------------------------------------------
function exportTSeriesMenuItem_Callback(hObject, eventdata, handles)
mrWarnDlg('exportTSeries not yet implemented');

% --------------------------------------------------------------------
function exportOverlayMenuItem_Callback(hObject, eventdata, handles)
pathstr = putPathStrDialog(pwd,'Specify name of Nifti file to export overlay to',mrGetPref('niftiFileExtension'));
% pathstr = [] if aborted
if ~isempty(pathstr)
  mrGlobals;
  viewNum = handles.viewNum;
  mrExport2SR(viewNum, pathstr);
end

% --------------------------------------------------------------------
function exportROIMenuItem_Callback(hObject, eventdata, handles)
% get view
mrGlobals;
viewNum = handles.viewNum;
v = viewGet(viewNum,'view');

% get the roi we are being asked to export
roiNum = viewGet(v,'currentroi');
if isempty(roiNum)
  mrWarnDlg('(mlrExportROI) No current ROI to export');
  return
end

% get current roi name
roiName = viewGet(v,'roiname');

% put up dialog to select filename
pathstr = putPathStrDialog(pwd,'Specify name of Nifti file to export ROI to',setext(roiName,mrGetPref('niftiFileExtension')));

% pathstr = [] if aborted
if ~isempty(pathstr)
  mlrExportROI(v, pathstr);
end


% --------------------------------------------------------------------
function exportImageMenuItem_Callback(hObject, eventdata, handles)
pathstr = putPathStrDialog(pwd,'Specify name of exported image file','*.tif');
% pathstr = [] if aborted
if ~isempty(pathstr)
    img = refreshMLRDisplay(handles.viewNum);
    imwrite(img, pathstr, 'tif');
end


% --------------------------------------------------------------------
function readmeMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
createReadme(MLR.session,MLR.groups);

% --------------------------------------------------------------------
function saveSessionMenuItem_Callback(hObject, eventdata, handles)
saveSession(1);

% --------------------------------------------------------------------
function printMenuItem_Callback(hObject, eventdata, handles)

mrGlobals;
viewNum = handles.viewNum;
v = MLR.views{viewNum};

mrPrint(v);

% --------------------------------------------------------------------
function quitMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
v = MLR.views{viewNum};

mrQuit(1,v);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function editMenu_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function editSessionMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
session = editSessionGUI('session',MLR.session);
% session is empty if aborted
if ~isempty(session)
    MLR.session = session;
    saveSession(0);
end

% --------------------------------------------------------------------
function editGroupMenu_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function newGroupMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
userInput = inputdlg('Enter name for new group: ','New group');
if ~isempty(userInput)
    groupName = userInput{1};
    view = viewSet(view,'newGroup',groupName);
end

% --------------------------------------------------------------------
function deleteGroupMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
groupNum = viewGet(view,'currentGroup');
view = viewSet(view,'deleteGroup',groupNum);
refreshMLRDisplay(viewNum);

% --------------------------------------------------------------------
function editGroupMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
curGroup = viewGet(view,'currentGroup');
group = viewGet(view,'group',curGroup);
if ~isempty(group.scanParams)
    group = editGroupGUI('group',MLR.groups(curGroup));
end
% group is empty if aborted
if ~isempty(group)
    MLR.groups(curGroup) = group;
    saveSession(0);
end

% --------------------------------------------------------------------
function infoGroupMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
groupNum = viewGet(view,'curGroup');
groupName = viewGet(view,'groupName',groupNum);
disp(sprintf('\n===== Group Info (%s) =====',groupName));
groupInfo(groupNum,0,1);
disp(sprintf('======================'));

% --------------------------------------------------------------------
function editScanMenu_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function addScanMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};

% Choose nifti file to add
tseriesDir = viewGet(view,'tseriesDir');
pathStr = mlrGetPathStrDialog(tseriesDir,'Add scan: choose nifti file',['*',mrGetPref('niftiFileExtension')]);
if isempty(pathStr)
    % Aborted
    return
end
if ~exist(pathStr,'file')
    mrErrorDlg('Invalid nifti file');
end

% Copy the nifti file to the tseries directory
[dir,file,ext] = fileparts(pathStr);
if strcmp(dir,tseriesDir)
    % If tseries file is already in the tseries directory, then use it
    hdr = mlrImageReadNiftiHeader(pathStr);
    fileName = [file,ext];
else
    % Copy file to the tseries directory
    fileName = ['tseries-',datestr(now,'mmddyy-HHMMSS'),mrGetPref('niftiFileExtension')];
    [data,hdr] = mlrImageReadNifti(pathStr);
    newPathStr = fullfile(tseriesDir,fileName);
    [bytes,hdr] = cbiWriteNifti(newPathStr,data,hdr);
end

% Add it
scanParams.fileName = fileName;
view = viewSet(view,'newScan',scanParams);

% Open GUI to set description, junk frames, nframes
curGroup = viewGet(view,'currentGroup');
group = editGroupGUI('group',MLR.groups(curGroup));
% group is empty if aborted
if ~isempty(group)
    MLR.groups(curGroup) = group;
    saveSession(0);
end

% --------------------------------------------------------------------
function copyScanMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
scanNum = viewGet(view,'currentScan');
scan.scanParams = viewGet(view,'scanParams',scanNum);
scan.scanParams.fileName = [];
% just save the filename instead of loading the whole tSeries
% since some scans are very large
%scan.tseries = loadTSeries(view,scanNum,'all');
scan.tseries = viewGet(view,'tseriesPathStr',scanNum);
MLR.clipboard = scan;

% --------------------------------------------------------------------
function pasteScanMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
if isfield(MLR.clipboard,'tseries') & isfield(MLR.clipboard,'scanParams') & ...
        isscan(MLR.clipboard.scanParams)
    view = saveNewTSeries(view,MLR.clipboard.tseries,MLR.clipboard.scanParams,MLR.clipboard.scanParams.niftiHdr);
else
    mrErrorDlg('(paste scan) Cannot paste. Clipboard does not contain a valid scan. Use Edit -> Scan -> Copy Scan.')
end

% --------------------------------------------------------------------
function deleteScanMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
if viewGet(view,'nScans') == 0
  disp(sprintf('(mrLoadRetGUI) No scans in group %s to delete',viewGet(view,'groupName')));
  return
end
deleteScans(view);

% --------------------------------------------------------------------
function transformsMenu_Callback(hObject, eventdata, handles)
%mrWarnDlg('transforms not yet implemented');

% --------------------------------------------------------------------
function sformScanMenuItem_Callback(hObject, eventdata, handles)

mrGlobals;
% get the view
viewNum = handles.viewNum;
v = MLR.views{viewNum};
% get the scanSform
scanSform = viewGet(v,'scansform');
sformCode = viewGet(v,'sformCode');
% params dialog
paramsInfo = {{'sform',scanSform,'The sform is usually set by mrAlign to specify the transformation from the scan coordinates to the volume anatomy coordinates. Only change this here if you know what you are doing! Also, any fix made here only changes the mrSession it does not change the original nifti header, so if you run mrUpdateNiftiHdr your change here will be overwritten.'},...
    {'sformCode',sformCode,'incdec=[-1 1]','minmax=[0 inf]','This gets set to 1 if mrAlign changes the sform. If it is 0 it means the sform has never been set. If you set this to 0 then mrLoadRet will ignore the sform as if it has never been set. If you want to change the above sform, make sure that this is 1'}};

params = mrParamsDialog(paramsInfo,'scanSform');

% ask the user if they are really sure before actually changing it
if ~isempty(params)
  answer = questdlg('Are you sure you want to change the sform (Normally you should fix problems with the sform by rerunning mrAlign. Also, any changes made here are only made to the mrSession variable they are not saved in the nifti header and will be overwritten if you ever call mrUpdateNifitHdr)?');
  if strcmp(answer,'Yes')
    v = viewSet(v,'scanSform',params.sform);
    v = viewSet(v,'sformCode',params.sformCode);
    saveSession;
    % clear the caches
    v = viewSet(v,'roiCache','init');
    v = viewSet(v,'overlayCache','init');
    v = viewSet(v,'baseCache','init');
    % and refresh
    refreshMLRDisplay(viewNum);
  end
end

% --------------------------------------------------------------------
function scan2BaseMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
% get the view
viewNum = handles.viewNum;
v = MLR.views{viewNum};
% get the scan2base
scan2base = inv(viewGet(v,'base2scan'));
sformCode = viewGet(v,'sformCode');
% params dialog
paramsInfo = {{'scan2base',scan2base,'This tells you the transformation from the scan coordinates to the base coordinates. If you have set the sfroms properly with mrAlign this should give an easily interpretable value. For instance if you have the same slices, but voxels are twice as big in the scan, then the diagonal elements should have 0.5 in them. This can be fixed here, but should only be done if you really know what you are doing. Otherwise this should be fixed by rerunning mrAlign. Also, any fix made here only changes the mrSession it does not change the original nifti header, so if you run mrUpdateNiftiHdr your change here will be overwritten.'},...
    {'sformCode',sformCode,'editable=0','This gets set to 1 if mrAlign changes the sform (3 for talairach). If it is 0 it means the sform has never been set. If you change the base2scan, it will automatically be set to 1.'}};

params = mrParamsDialog(paramsInfo,'scan2base transformation');

if ~isempty(params)
  answer = questdlg('Are you sure you want to change the sform (Normally you should fix problems with the sform by rerunning mrAlign. Also, any changes made here are only made to the mrSession variable they are not saved in the nifti header and will be overwritten if you ever call mrUpdateNifitHdr)?');
  if strcmp(answer,'Yes')
    v = viewSet(v,'scan2base',params.scan2base);
    saveSession;
    % clear the caches
    v = viewSet(v,'roiCache','init');
    v = viewSet(v,'overlayCache','init');
    v = viewSet(v,'baseCache','init');
    % and refresh
    refreshMLRDisplay(viewNum);
  end
end

% --------------------------------------------------------------------
function dicomInfoMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
% get the view
viewNum = handles.viewNum;
v = MLR.views{viewNum};
s = viewGet(v,'curScan');
g = viewGet(v,'curGroup');
dicomInfo(s,g);

% --------------------------------------------------------------------
function infoScanMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
scanNum = viewGet(view,'curScan');
groupNum = viewGet(view,'curGroup');
groupName = viewGet(view,'groupName',groupNum);
scanInfo(scanNum,groupNum,1);


% --------------------------------------------------------------------
function linkStimfileMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
% Choose stimfile to add
etcDir = viewGet(view,'etcDir');
pathStr = mlrGetPathStrDialog(etcDir,'Link Stimfile: choose matlab file');
if isempty(pathStr)
    % Aborted
    return
end
if ~exist(pathStr,'file')
    mrErrorDlg('Invalid stimfile file');
end

% Copy the nifti file to the tseries directory
[dir,file,ext] = fileparts(pathStr);
if strcmp(dir,etcDir)
    % If stimfile is already in the Etc directory
    fileName = [file,ext];
    viewSet(view, 'stimfilename', fileName);
else
    mrErrorDlg('Stimfile must be in the Etc directory');
end

% --------------------------------------------------------------------
function editAnalysisMenu_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function newAnalysisMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
userInput = inputdlg('Enter name for new analysis: ','New analysis');
if ~isempty(userInput)
  newAnalysis(view,userInput{1});
  refreshMLRDisplay(viewNum);
end

% --------------------------------------------------------------------
function copyAnalysisMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
MLR.clipboard = viewGet(view,'analysis');

% --------------------------------------------------------------------
function pasteAnalysisMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
[check analysis] = isanalysis(MLR.clipboard);
if check
    view = viewSet(view,'newAnalysis',analysis);
else
    mrErrorDlg('(paste analysis) Cannot paste. Clipboard does not contain a valid analysis. Use Edit -> Analysis -> Copy Analysis.')
end
refreshMLRDisplay(viewNum);

% --------------------------------------------------------------------
function editAnalysisMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
%view = editAnalysisGUI(view);

paramsInfo{1} = {'name',viewGet(view,'analysisName'),'Name of analysis'};
paramsInfo{end+1} = {'date',viewGet(view,'analysisDate'),'editable=0','Date analysis was done (not editable)'};
paramsInfo{end+1} = {'function',viewGet(view,'analysisFunction'),'editable=0','Function with which analysis was done (not editable)'};
paramsInfo{end+1} = {'groupName',viewGet(view,'analysisGroupName'),'editable=0','Group on which analysis was done (not editable)'};
paramsInfo{end+1} = {'GUIFunction',viewGet(view,'analysisGUIFunction'),'GUI function for analysis'};
paramsInfo{end+1} = {'reconcileFunction',viewGet(view,'analysisReconcileFunction'),'Reconcile function for analysis'};
paramsInfo{end+1} = {'mergeFunction',viewGet(view,'analysisMergeFunction'),'Merge function for analysis'};
paramsInfo{end+1} = {'type',viewGet(view,'analysisType'),'Type of analysis'};
paramsInfo{end+1} = {'clipAcrossOverlays',viewGet(view,'clipAcrossOverlays'),'type=checkbox','Clipping behaviour: if 1, only voxels that are inside clipping values in all overlays are displayed; if 0, only the clipping values of the current overlay and its alpha overlay are taken into account'};

params = mrParamsDialog(paramsInfo,'Edit analysis');
if ~isempty(params)
  view = viewSet(view,'analysisName',params.name);
  view = viewSet(view,'analysisGUIFunction',params.GUIFunction);
  view = viewSet(view,'analysisReconcileFunction',params.reconcileFunction);
  view = viewSet(view,'analysisMergeFunction',params.mergeFunction);
  view = viewSet(view,'analysisType',params.type);
  view = viewSet(view,'clipAcrossOverlays',params.clipAcrossOverlays);
end
refreshMLRDisplay(viewNum);

% --------------------------------------------------------------------
function EditAnalysisInfoMenuItem_Callback(hObject, eventdata, handles)

mrGlobals;
viewNum = handles.viewNum;
v = MLR.views{viewNum};

% no current anatomy, just return
if isempty(viewGet(v,'curAnalysis')),return;end

disppercent(-inf,'Gathering analysis info');
% get the current analysis
a = viewGet(v,'Analysis',viewGet(v,'curAnalysis'));

% get fields
fields = fieldnames(a);
fields = setdiff(fields,{'d','overlays','params','curOverlay'});

% make into a display
paramsInfo = {};
for fnum = 1:length(fields)
  if ~isstruct(fields)
    paramsInfo{end+1} = {fields{fnum},a.(fields{fnum}),'editable=0'};
  end
end

% check d
if isfield(a,'d')
  dExists = [];
  for dnum = 1:length(a.d)
    dExists(dnum) = ~isempty(a.d{dnum});
  end
  paramsInfo{end+1} = {sprintf('dScans'),num2str(find(dExists)),'editable=0',sprintf('Scans that d structure exists for')};
end

% check overlays
for onum = 1:length(a.overlays)
  for snum = 1:length(a.overlays(onum).data)
    overlayExists(snum) = ~isempty(a.overlays(onum).data{snum});
  end
  % make params for this
  paramsInfo{end+1} = {sprintf('overlay%i',onum),a.overlays(onum).name,'editable=0',sprintf('Name of overlay %i',onum)};
  paramsInfo{end+1} = {sprintf('overlay%iScans',onum),num2str(find(overlayExists)),'editable=0',sprintf('Scans that overlay %s exists for',a.overlays(onum).name)};
end

if isfield(a.params,'paramInfo') %if the parameteres have been asked by mrParams
  paramsInfo{end+1} = {'params',[],'View analysis parameters','type=pushbutton','buttonString=View analysis parameters','callback',{@mrParamsDisp,a.params,'Analysis Parameters'}};
else
  paramsInfo{end+1} = {'params',[],'View analysis parameters','type=pushbutton','buttonString=View analysis parameters','callback',@viewAnalysisParams,'callbackArg',v};
end

disppercent(inf);

% display parameters
mrParamsDialog(paramsInfo,'Analysis Info');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% helper function to view params
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function retval = viewAnalysisParams(v)

% bogus return value
retval = [];

% get analysis info
curAnalysis = viewGet(v,'curAnalysis');
params = viewGet(v,'analysisParams',curAnalysis);
guiFunction = viewGet(v,'analysisGuiFunction',curAnalysis);
groupName = viewGet(v,'analysisGroupName',curAnalysis);

% check for function
while exist(sprintf('%s.m',stripext(guiFunction))) ~= 2
  paramsInfo = {{'GUIFunction',guiFunction,sprintf('The GUI function for this analysis, %s, was not found. If you want to specify another function name you can enter that here and try again.',guiFunction)}};
  paramsGUIFunction = mrParamsDialog(paramsInfo,sprintf('GUI function: %s not found',guiFunction));
  if isempty(paramsGUIFunction) || strcmp(paramsGUIFunction.GUIFunction,guiFunction)
    return
  else
    guiFunction = paramsGUIFunction.GUIFunction;
  end
end

% params = guiFunction('groupName',groupName,'params',params);
evalstring = ['params = ',guiFunction,'(','''','groupName','''',',groupName,','''','params','''',',params);'];

eval(evalstring);


% --------------------------------------------------------------------
function editOverlayMenu_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function copyOverlayMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
MLR.clipboard = viewGet(view,'overlay');

% --------------------------------------------------------------------
function pasteOverlayMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
[check overlay] = isoverlay(MLR.clipboard);
if ~check
    mrErrorDlg('(paste overlay) Cannot paste. Clipboard does not contain a valid overlay. Use Edit -> Overlay -> Copy Overlay.')
end
if ~isanalysis(viewGet(view,'analysis'))
    mrErrorDlg('(paste overlay) Overlays must be pasted into an analysis. Use Edit -> Analysis -> New Analysis.')
end
view = viewSet(view,'newOverlay',overlay);
refreshMLRDisplay(viewNum);

% --------------------------------------------------------------------
function editOverlayMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
editOverlayGUImrParams(viewNum);
% view = MLR.views{viewNum};
% view = editOverlayGUI(view);
% view = viewSet(view,'overlayCache','init');
% refreshMLRDisplay(viewNum);

% --------------------------------------------------------------------
function overlayInfoMenuItem_Callback(hObject, eventdata, handles)

mrGlobals;
viewNum = handles.viewNum;
v = MLR.views{viewNum};

overlayInfo(v);

% --------------------------------------------------------------------
function editRoiMenu_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function copyRoiMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
MLR.clipboard = viewGet(view,'roi');

% --------------------------------------------------------------------
function pasteRoiMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
refresh=false;
for iRoi = 1:length(MLR.clipboard)
  % Check to see that it is a valid ROI structure and then add it.
  [check roi] = isroi(MLR.clipboard(iRoi));
  if check
      view = viewSet(view,'newROI',roi);
      refresh=true;
  else
      mrWarnDlg('(paste ROI) Cannot paste. Clipboard does not contain a valid ROI. Use Edit -> ROI -> Copy ROI.')
  end
end
if refresh
  % Select the last ROI
  ROInum = viewGet(view,'numberofROIs');
  if (ROInum > 0)
      view = viewSet(view,'currentROI',ROInum);
  end
  refreshMLRDisplay(viewNum);
end

% --------------------------------------------------------------------
function editRoiMenuItem_Callback(hObject, eventdata, handles)

viewNum = handles.viewNum;
editManyROIs(viewNum,viewGet(viewNum,'currentROI'));

% --------------------------------------------------------------------
function editManyROIsMenuItem_Callback(hObject, eventdata, handles)
viewNum = handles.viewNum;
editManyROIs(viewNum,1:viewGet(viewNum,'numROIs'));

function editManyROIs(viewNum,roiList)

if isempty(roiList)
  disp('(editManyROIs) No ROI is loaded')
  return;
end
mrGlobals;
v = MLR.views{viewNum};
paramsInfo = {};
for roinum = roiList
  % get name and colors for each roi
  roiNames{roinum} = viewGet(v,'roiName',roinum);
  roiNotes = viewGet(v,'roiNotes',roinum);
  colors = putOnTopOfList(viewGet(v,'roiColor',roinum),color2RGB);
  displayOn = putOnTopOfList(viewGet(v,'roiDisplayOnBase'),viewGet(v,'baseNames'));
  paramsInfo{end+1} = {sprintf('%sName',fixBadChars(roiNames{roinum})),roiNames{roinum},'Name of roi, avoid using punctuation and space'};
  paramsInfo{end+1} = {sprintf('%sColor',fixBadChars(roiNames{roinum})),colors,'type=popupmenu',sprintf('The color that roi %s will display in',roiNames{roinum})};
  paramsInfo{end+1} = {sprintf('%sNotes',fixBadChars(roiNames{roinum})),roiNotes,sprintf('Note for roi %s',roiNames{roinum})};
  paramsInfo{end+1} = {sprintf('%sDisplayOnBase',fixBadChars(roiNames{roinum})),displayOn,sprintf('Base that roi %s is best displayed on',roiNames{roinum})};
end
if isempty(paramsInfo),return,end
params = mrParamsDialog(paramsInfo,'Edit Many ROIs');

% if not empty, then change the parameters
if ~isempty(params)
  for roinum = roiList
    roiName = fixBadChars(roiNames{roinum});
    v = viewSet(v,'roiColor',params.(sprintf('%sColor',roiName)),roinum);
    v = viewSet(v,'roiName',params.(sprintf('%sName',roiName)),roinum);
    v = viewSet(v,'roiNotes',params.(sprintf('%sNotes',roiName)),roinum);
    v = viewSet(v,'roiDisplayOnBase',params.(sprintf('%sDisplayOnBase',roiName)),roinum);
  end
  refreshMLRDisplay(viewNum);
end



% --------------------------------------------------------------------
function editAllROIsMenuItem_Callback(hObject, eventdata, handles)

mrGlobals;
viewNum = handles.viewNum;
v = MLR.views{viewNum};
nROIs = viewGet(v,'numROIs');

if nROIs == 0,return,end
paramsInfo = {};
colors = putOnTopOfList(viewGet(v,'roiColor'),color2RGB);
% get color to set all ROIs to
paramsInfo{end+1} = {'changeColor',0,'type=checkbox','Click to change color of all ROIs'};
paramsInfo{end+1} = {'color',colors,'type=popupmenu','contingent=changeColor','Color for all ROIs'};
paramsInfo{end+1} = {'changeNotes',0,'type=checkbox','Click to change notes of all ROIs'};
paramsInfo{end+1} = {'notes',viewGet(v,'roiNotes'),'contingent=changeNotes','Notes for all ROIs'};
params = mrParamsDialog(paramsInfo,'Edit All ROIs',1.5);

% if not empty, then change the parameters
if ~isempty(params)
  for roinum = 1:nROIs
    if ~isempty(params.color)
      v = viewSet(v,'roiColor',params.color,roinum);
    end
    if ~isempty(params.notes)
      v = viewSet(v,'roiNotes',params.notes,roinum);
    end
  end
  refreshMLRDisplay(viewNum);
end


% --------------------------------------------------------------------
function infoROIMenuItem_Callback(hObject, eventdata, handles)

mrGlobals;
viewNum = handles.viewNum;
v = MLR.views{viewNum};

roiInfo(v);

% --------------------------------------------------------------------
function editBaseMenu_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function copyBaseMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
MLR.clipboard = viewGet(view,'baseAnatomy');

% --------------------------------------------------------------------
function pasteBaseMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
[check base] = isbase(MLR.clipboard);
if check
    view = viewSet(view,'newBase',base);
else
    mrErrorDlg('(paste base anatomy) Cannot paste. Clipboard does not contain a valid scan. Use Edit -> Base Anatomy -> Copy Base Anatomy.')
end
refreshMLRDisplay(viewNum);

% --------------------------------------------------------------------
function editBaseMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
v = MLR.views{viewNum};
v = editBaseGUI(v);
refreshMLRDisplay(viewNum);

% --------------------------------------------------------------------
function baseTransformsMenu_Callback(hObject, eventdata, handles)
% hObject    handle to baseTransformsMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function baseSformMenuItem_Callback(hObject, eventdata, handles)

mrGlobals;
% get the view
viewNum = handles.viewNum;
v = MLR.views{viewNum};
% get the baseSform
baseSform = viewGet(v,'baseSform');
sformCode = viewGet(v,'basesformCode');
% params dialog
paramsInfo = {{'sform',baseSform,'The sform is usually set by mrAlign to specify the transformation from the base coordinates to the volume anatomy coordinates. Only change this here if you really know what you are doing!'},...
    {'sformCode',sformCode,'incdec=[-1 1]','minmax=[0 inf]','This gets set to 1 if mrAlign changes the sform. If it is 0 it means the sform has never been set. If you set this to 0 then mrLoadRet will ignore the sform as if it has never been set. If you want to change the above sform, make sure that this is 1'}};

params = mrParamsDialog(paramsInfo,'baseSform');

% ask the user if they are really sure before actually changing it
if ~isempty(params)
  answer = questdlg('Are you sure you want to change the sform (Normally you should fix problems with the sform by rerunning mrAlign.');
  if strcmp(answer,'Yes')
    v = viewSet(v,'baseSform',params.sform);
    v = viewSet(v,'baseSformCode',params.sformCode);
    saveSession;
    % clear the caches
    v = viewSet(v,'roiCache','init');
    v = viewSet(v,'overlayCache','init');
    v = viewSet(v,'baseCache','init');
    % and refresh
    refreshMLRDisplay(viewNum);
  end
end


% --------------------------------------------------------------------
function base2scanMenuItem_Callback(hObject, eventdata, handles)

mrGlobals;
% get the view
viewNum = handles.viewNum;
v = MLR.views{viewNum};
% get the base2scan xform
base2scan = viewGet(v,'base2scan');
baseSformCode = viewGet(v,'baseSformCode');
% params dialog
paramsInfo = {{'base2scan',base2scan,'This tells you the transformation from the base coordinates to the scan coordinates. If you have set the sfroms properly with mrAlign this should give an easily interpretable value. For instance if you have the same slices, but voxels are twice as big in the scan, then the diagonal elements should have 0.5 in them. This can be fixed here, but should only be done if you really know what you are doing. Otherwise this should be fixed by rerunning mrAlign.'},...
    {'baseSformCode',baseSformCode,'editable=0','This gets set to 1 if mrAlign changes the sform (3 for talairach). If it is 0 it means the sform has never been set. If you change the base2scan, it will automatically be set to 1.'}};

params = mrParamsDialog(paramsInfo,'base2scan transformation');

if ~isempty(params)
  answer = questdlg('Are you sure you want to change the sform (Normally you should fix problems with the sform by rerunning mrAlign.');
  if strcmp(answer,'Yes')
    v = viewSet(v,'base2scan',params.base2scan);
    saveSession;
    % clear the caches
    v = viewSet(v,'roiCache','init');
    v = viewSet(v,'overlayCache','init');
    v = viewSet(v,'baseCache','init');
    % and refresh
    refreshMLRDisplay(viewNum);
  end
end


% --------------------------------------------------------------------
function infoBaseAnatomyMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
baseInfo(view);

% --------------------------------------------------------------------
function prefMenu_Callback(hObject, eventdata, handles)

mrGlobals;
viewNum = handles.viewNum;
v = MLR.views{viewNum};

% get settings from start
interpMethod = mrGetPref('interpMethod');
selectedROIColor = mrGetPref('selectedROIColor');
roiContourWidth = mrGetPref('roiContourWidth');
corticalDepthBins = mrGetPref('corticalDepthBins');
roiCorticalDepthDisplayRatio = mrGetPref('roiCorticalDepthDisplayRatio');

% remember old cache sizes
roiCacheSize = mrGetPref('roiCacheSize');
baseCacheSize = mrGetPref('baseCacheSize');
overlayCacheSize = mrGetPref('overlayCacheSize');

% prefs dialog
prefParams = mrEditPrefs;

% if changes, check for change in cache size
if ~isempty(prefParams)
  if (roiCacheSize ~= prefParams.roiCacheSize)
    v = viewSet(v,'roiCache','init');
  end
  if (baseCacheSize ~= prefParams.baseCacheSize)
    v = viewSet(v,'baseCache','init');
  end
  if (overlayCacheSize ~= prefParams.overlayCacheSize)
    v = viewSet(v,'overlayCache','init');
  end

  % check to see if interpolation method has changed
  if ~strcmp(interpMethod,prefParams.interpMethod)
    % dump overlay cache and redraw
    v = viewSet(v,'overlayCache','init');
    refreshMLRDisplay(v.viewNum);
  end
  
  % check to see if parameters affecting the ROI cache have changed
  if ~isequal(roiCorticalDepthDisplayRatio,prefParams.roiCorticalDepthDisplayRatio)
    % dump overlay cache and redraw
    v = viewSet(v,'roiCache','init');
    refreshMLRDisplay(v.viewNum);
  end
  
  if ~strcmp(corticalDepthBins,prefParams.corticalDepthBins)
    set(handles.corticalDepthSlider,'SliderStep',min(1/(prefParams.corticalDepthBins-1)*[1 3],1));
    if isfield(handles,'corticalMaxDepthSlider')
      set(handles.corticalMaxDepthSlider,'SliderStep',min(1/(prefParams.corticalDepthBins-1)*[1 3],1));
    end
    refreshMLRDisplay(v.viewNum);
  end
  
  % check to see if any other parameters that affects the display, but not the cache have changed
  if ~strcmp(selectedROIColor,prefParams.selectedROIColor) ||...
      ~strcmp(roiContourWidth,prefParams.roiContourWidth) || ...
    refreshMLRDisplay(v.viewNum);
  end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function windowMenu_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function graphMenuItem_Callback(hObject, eventdata, handles)
newGraphWin;

% --------------------------------------------------------------------
function newWindowMenuItem_Callback(hObject, eventdata, handles)
view = mrOpenWindow('',[]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function analysisMenu_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function motionCompMenu_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function motionCompMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
view = motionComp(view);

% --------------------------------------------------------------------
function motionCompWithinMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
view = motionCompWithinScan(view);

% --------------------------------------------------------------------
function motionCompBetweenMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
view = motionCompBetweenScans(view);

% --------------------------------------------------------------------
function sliceTimeCorrectMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
view = sliceTimeCorrect(view);

% --------------------------------------------------------------------
function averageTSeriesMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
view = averageTSeries(view);

% --------------------------------------------------------------------
function concatenateTSeriesMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
view = concatTSeries(view);

% --------------------------------------------------------------------
function tsStatsMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
view = timeSeriesStats(view);
refreshMLRDisplay(viewNum);

% --------------------------------------------------------------------
function corAnalMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
view = corAnal(view);
refreshMLRDisplay(viewNum);

% --------------------------------------------------------------------
function eventRelatedMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
view = eventRelated(view);

% --------------------------------------------------------------------
function glmMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
view = eventRelatedGlm(view);

% --------------------------------------------------------------------
function recomputeAnalysisMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
if viewGet(view,'numAnalyses') > 0
  n = viewGet(view,'currentAnalysis');
  groupName = viewGet(view,'analysisGroupName',n);
  analysisFunction = viewGet(view,'analysisFunction',n);
  guiFunction = viewGet(view,'analysisGuiFunction',n);
  params = viewGet(view,'analysisParams',n);
  % run parameters again if guiFunction is valid
  if ~isempty(guiFunction)
    % params = guiFunction('groupName',groupName,'params',params);
    evalstring = ['params = ',guiFunction,'(','''','groupName','''',',groupName,','''','params','''',',params);'];
    eval(evalstring);
  else
    mrWarnDlg('(mrLoadRetGUI) No guiFunction set for analysis to display parameters');
  end
  % params is empty if GUI cancelled
  if ~isempty(params) && ~isempty(analysisFunction)
    % view = analysisFunction(view,params);
    evalstring = ['view = ',analysisFunction,'(view,params);'];
    eval(evalstring);
    refreshMLRDisplay(viewNum);
  elseif isempty(analysisFunction)
    mrWarnDlg('(mrLoadRetGUI) No analysis function set for analysis to rerun analysis');
  end

else
  mrWarnDlg(sprintf('(mrLoadRetGUI) No analyses loaded'));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function viewMenu_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function deleteBaseMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
baseNum = viewGet(view,'currentBase');
if ~isempty(baseNum)
  view = viewSet(view,'deleteBase',baseNum);
  refreshMLRDisplay(viewNum);
end

% --------------------------------------------------------------------
function deleteAllBasesMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
numBases = viewGet(view,'numberofbasevolumes');
for baseNum = numBases:-1:1;
    view = viewSet(view,'deleteBase',baseNum);
end
refreshMLRDisplay(viewNum);

% --------------------------------------------------------------------
function deleteManyBasesMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
numBases =  selectInList(view,'bases','Select bases to remove');
if ~isempty(numBases)
  for baseNum = fliplr(numBases);
      view = viewSet(view,'deleteBase',baseNum);
  end
  refreshMLRDisplay(viewNum);
end

% --------------------------------------------------------------------
function deleteAnalysisMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
analysisNum = viewGet(view,'currentAnalysis');
view = viewSet(view,'deleteAnalysis',analysisNum);
refreshMLRDisplay(viewNum);

% --------------------------------------------------------------------
function deleteManyAnalysisMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
numAnalyses =  selectInList(view,'analyses','Select analyses to remove');
if ~isempty(numAnalyses)
   view = viewSet(view,'deleteAnalysis',numAnalyses);
   refreshMLRDisplay(viewNum);
end

% --------------------------------------------------------------------
function deleteAllAnalysisMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
numAnalyses = viewGet(view,'numberofAnalyses');
view = viewSet(view,'deleteAnalysis',1:numAnalyses);
% also, go through and delete all "laodedAnalyses", i.e.
% analysis that are loaded in ohter groups
for g = 1:viewGet(view,'numGroups');
  view = viewSet(view,'loadedAnalyses',[],g);
end
refreshMLRDisplay(viewNum);

% --------------------------------------------------------------------
function deleteOverlayMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
overlayNum = viewGet(view,'currentOverlay');
view = viewSet(view,'deleteOverlay',overlayNum);
refreshMLRDisplay(viewNum);

% --------------------------------------------------------------------
function deleteManyOverlaysMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
numOverlays = selectInList(view,'overlays','Select overlays to remove');
if ~isempty(numOverlays)
   view = viewSet(view,'deleteOverlay',numOverlays);
   refreshMLRDisplay(viewNum);
end

% --------------------------------------------------------------------
function deleteAllOverlayMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
numOverlays = viewGet(view,'numberofoverlays');
view = viewSet(view,'deleteOverlay',1:numOverlays);
refreshMLRDisplay(viewNum);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function roiMenu_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function createRoiMenu_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function newRoiMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
[view userCancel] = newROI(view);
if userCancel,return,end
refreshMLRDisplay(viewNum);

% --------------------------------------------------------------------
function createRectangleMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
[view userCancel] = newROI(view);
if userCancel,return,end
view = drawROI(view,'rectangle',1);
refreshMLRDisplay(viewNum);

% --------------------------------------------------------------------
function createPolygonMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
[view userCancel] = newROI(view);
if userCancel,return,end
view = drawROI(view,'polygon',1);
refreshMLRDisplay(viewNum);

% --------------------------------------------------------------------
function createLineMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
[view userCancel] = newROI(view);
if userCancel,return,end
view = drawROI(view,'line',1);
refreshMLRDisplay(viewNum);

% --------------------------------------------------------------------
function createContiguousMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
[view userCancel] = newROI(view);
if userCancel,return,end
view = drawROI(view,'contiguous',1);
refreshMLRDisplay(viewNum);

% --------------------------------------------------------------------
function deleteRoiMenu_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function deleteRoiMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
roinum = viewGet(view,'currentROI');
view = viewSet(view,'deleteROI',roinum);
refreshMLRDisplay(viewNum);

% --------------------------------------------------------------------
function removeManyROIMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
thisView = MLR.views{viewNum};

% put up a dialog with rois to delete
roiNums = selectInList(thisView,'rois','Select ROIs to remove');
if ~isempty(roiNums)
  % now delete anything the user selected
  thisView = viewSet(thisView,'deleteROI',roiNums);
  refreshMLRDisplay(viewNum);
end

% --------------------------------------------------------------------
function deleteAllROIsMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
numrois = viewGet(view,'numberofrois');
for roinum = numrois:-1:1;
    view = viewSet(view,'deleteROI',roinum);
end
refreshMLRDisplay(viewNum);

% --------------------------------------------------------------------
function addRoiMenu_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function addRectangleMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
view = drawROI(view,'rectangle',1);
refreshMLRDisplay(viewNum);

% --------------------------------------------------------------------
function addPolygonMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
view = drawROI(view,'polygon',1);
refreshMLRDisplay(viewNum);

% --------------------------------------------------------------------
function addLineMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
view = drawROI(view,'line',1);
refreshMLRDisplay(viewNum);


% --------------------------------------------------------------------
function addContiguousMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
view = drawROI(view,'contiguous',1);
refreshMLRDisplay(viewNum);


% --------------------------------------------------------------------
function removeRoiMenu_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function removeRectangleMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
view = drawROI(view,'rectangle',0);
refreshMLRDisplay(viewNum);

% --------------------------------------------------------------------
function removePolygonMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
view = drawROI(view,'polygon',0);
refreshMLRDisplay(viewNum);

% --------------------------------------------------------------------
function removeLineMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
view = drawROI(view,'line',0);
refreshMLRDisplay(viewNum);

% --------------------------------------------------------------------
function removeContiguousMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
view = drawROI(view,'contiguous',0);
refreshMLRDisplay(viewNum);

% --------------------------------------------------------------------
function combineROIMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
roiNames = viewGet(view,'roiNames');
selectedRoiNames = viewGet(view,'roiName');
if ischar(selectedRoiNames) || isempty(selectedRoiNames)
  roiNames1 = putOnTopOfList(selectedRoiNames,roiNames);
else
  roiNames1 = putOnTopOfList(selectedRoiNames{1},roiNames);
end
if iscell(selectedRoiNames) && length(selectedRoiNames)>1
  roiNames2 = putOnTopOfList(selectedRoiNames{2},roiNames);
else
  roiNames2=roiNames;
end
paramInfo = {...
  {'combineROI',roiNames1,'type=popupmenu','The ROI that will be combined with the otherROI'},...
  {'otherROI',roiNames2,'type=popupmenu','The otherROI is combined with the combineROI and the result is put into combineROI.'},...
  {'action',{'A not B', 'Intersection', 'Union', 'XOR'},'type=popupmenu','Select action for combining ROIs.'},...
  {'newName','','Add a name here if you want to make the combination have a new roi name. Otherwise it will rewrite the combineROI.'},...
  {'combine',0,'type=pushbutton','callback',@doCombine,'passParams=1','callbackArg',viewNum,'buttonString=Do combination','Click this button to do the combination. This is the same as hitting OK but won''t close the dialog so you can continue to do more combinations'}};
params = mrParamsDialog(paramInfo,'Combine ROIs');
if ~isempty(params)
  doCombine(viewNum,params);
end

function retval = doCombine(viewNum,params)

% get the view (it is important to get it from the viewNum
% so that we always have an up-to-date view
v = viewGet([],'view',viewNum);

retval = 1;
if isempty(params.newName)
  disp(sprintf('(doCombine) %s %s %s',params.combineROI,params.action,params.otherROI));
  v = combineROIs(v,params.combineROI,params.otherROI,params.action);
else
  disp(sprintf('(doCombine) %s %s %s: Saving as %s',params.combineROI,params.action,params.otherROI,params.newName));
  v = combineROIs(v,params.combineROI,params.otherROI,params.action,params.newName);
end
refreshMLRDisplay(viewNum);

% --------------------------------------------------------------------
function restrictRoiMenu_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function restrictRoiMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
currentROIs = viewGet(view,'currentROI');
scan = viewGet(view,'curscan');
for roinum = currentROIs
    view = restrictROI(view,roinum,scan);
end
refreshMLRDisplay(viewNum);

% --------------------------------------------------------------------
function restrictAllROIsMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
scan = viewGet(view,'curscan');
nROIs = viewGet(view,'numberofROIs');
for roinum = 1:nROIs
    view = restrictROI(view,roinum,scan);
end
refreshMLRDisplay(viewNum);

% --------------------------------------------------------------------
function undoRoiMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
prevCoords = viewGet(view,'prevROIcoords');
curCoords = viewGet(view,'roiCoords');
if ischar(prevCoords) || isequal(prevCoords,curCoords)
  disp('(mrLoadRetGUI: undoROI) Nothing to undo');
  return
end
if length(viewGet(view,'currentRoi'))>1
  disp('(mrLoadRetGUI: undoROI) Undo not implemented for several selected ROIs');
  return;
end
view = viewSet(view,'prevROIcoords',curCoords);
view = viewSet(view,'ROIcoords',prevCoords);
refreshMLRDisplay(viewNum);


% --------------------------------------------------------------------
function convertCorticalDepthRoiMenuItem_Callback(hObject, eventdata, handles)

mrGlobals;
viewNum = handles.viewNum;
v = MLR.views{viewNum};

v = convertROICorticalDepth(v);


% --------------------------------------------------------------------
function convertRoiCoordsMenuItem_Callback(hObject, eventdata, handles)

mrGlobals;
viewNum = handles.viewNum;
v = MLR.views{viewNum};

v = convertROI(v);


% --------------------------------------------------------------------
function showRoiMenu_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function showAllMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
view = viewSet(view,'showROIs','all');
refreshMLRDisplay(viewNum);

% --------------------------------------------------------------------
function showAllPerimeterMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
view = viewSet(view,'showROIs','all perimeter');
refreshMLRDisplay(viewNum);

% --------------------------------------------------------------------
function showSelectedMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
view = viewSet(view,'showROIs','selected');
refreshMLRDisplay(viewNum);

% --------------------------------------------------------------------
function showSelectedPerimeterMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
view = viewSet(view,'showROIs','selected perimeter');
refreshMLRDisplay(viewNum);

% --------------------------------------------------------------------
function setROIGroupMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to setROIGroupMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

mrGlobals;
viewNum = handles.viewNum;
v = MLR.views{viewNum};

inGroup = [];
roiGroup = viewGet(v,'roiGroup');

% get number of ROIs
nROIs = viewGet(v,'nrois');

% get roi names
if nROIs > 0
  roiNames = viewGet(v,'roiNames');
  roiGroupNum = zeros(1,nROIs);
  if ~isempty(roiGroup)
    roiGroupNum(roiGroup) = 1;
  end
  for iROI = 1:nROIs
    paramsInfo{iROI} = {roiNames{iROI},roiGroupNum(iROI),'type=checkbox'};
  end
  paramsInfo{end+1} = {'all',0,'type=pushbutton','callback',@setSpecialGroupAll,'buttonString','all','passParams=1','callbackArg',{roiNames 1},'Set all of the ROIS to be shown in the group'};
  paramsInfo{end+1} = {'none',0,'type=pushbutton','callback',@setSpecialGroupAll,'buttonString','none','passParams=1','callbackArg',{roiNames 0},'Set none of the ROIs to be shown in the group'};
  paramsInfo{end+1} = {'begins','','type=string','callback',@setSpecialGroupBeginsWith,'passParams=1','callbackArg',roiNames,'Set all the ROIs that begin with the letters input'};
  paramsInfo{end+1} = {'contains','','type=string','callback',@setSpecialGroupContains,'passParams=1','callbackArg',roiNames,'Set all the ROIs that contain the letters input'};
  
  params = mrParamsDialog(paramsInfo,'Choose ROIs in group',0.5);
  if isempty(params),return,end
  
  % create a cell array with the names of the ROIs to display
  roiGroup = {};
  for i = 1:length(roiNames)
    if params.(fixBadChars(roiNames{i}))
      roiGroup{end+1} = roiNames{i};
    end
  end

  % set the roi group
  v = viewSet(v,'roiGroup',roiGroup);

  % if the setting is not to show groups, change it to show goups
  showROIs = viewGet(v,'showROIs');
  if ~any(strcmp(showROIs,{'group','group perimeter'}))
    % see if perimeter is selected
    if ~isempty(strfind(showROIs,'perimeter'))
      v = viewSet(v,'showROIs','group perimeter');
    else
      v = viewSet(v,'showROIs','group');
    end
  end
  refreshMLRDisplay(viewGet(v,'viewNum'));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    setSpecialGroupAll    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function retval = setSpecialGroupAll(arg,params)

retval = 0;

% set all of them
for iROI = 1:length(arg{1})
  params.(arg{1}{iROI}) = arg{2};
end

params = mrParamsSet(params);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    setSpecialGroupBeginsWith    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function retval = setSpecialGroupBeginsWith(roiNames,params)

retval = 0;

% set all of them
for iROI = 1:length(roiNames)
  if strncmp(roiNames{iROI},params.begins,length(params.begins))
    params.(roiNames{iROI}) = 1;
  else
    params.(roiNames{iROI}) = 0;
  end
end

params = mrParamsSet(params);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    setSpecialGroupContains    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function retval = setSpecialGroupContains(roiNames,params)

retval = 0;

% set all of them
for iROI = 1:length(roiNames)
  if ~isempty(strfind(roiNames{iROI},params.contains))
    params.(roiNames{iROI}) = 1;
  else
    params.(roiNames{iROI}) = 0;
  end
end

params = mrParamsSet(params);

% --------------------------------------------------------------------
function showGroupMenuItem_Callback(hObject, eventdata, handles)

mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
view = viewSet(view,'showROIs','group');
refreshMLRDisplay(viewNum);

% --------------------------------------------------------------------
function showGroupPerimeterMenuItem_Callback(hObject, eventdata, handles)

mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
view = viewSet(view,'showROIs','group perimeter');
refreshMLRDisplay(viewNum);
  
% --------------------------------------------------------------------
function hideROIsMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
view = viewSet(view,'showROIs','hide');
refreshMLRDisplay(viewNum);

% --------------------------------------------------------------------
function labelsROIsMenuItem_Callback(hObject, eventdata, handles)

mrGlobals;
viewNum = handles.viewNum;
v = MLR.views{viewNum};

% switch labels on/off
labelROIs = viewGet(v,'labelROIs');
viewSet(v,'labelROIs',~labelROIs');

% redisplay
refreshMLRDisplay(viewNum);

% --------------------------------------------------------------------
function findCurrentROIMenuItem_Callback(hObject, eventdata, handles)

mrGlobals;
viewNum = handles.viewNum;
v = MLR.views{viewNum};

v = findROI(v);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotMenu_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function interrogateOverlayMenuItem_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
% start or stop the interrogator
if strcmp(get(hObject,'Checked'),'on')
    mrInterrogator('end',viewNum);
    set(hObject,'Checked','off');
else
    % start the mrInterrogator
    mrInterrogator('init',viewNum);
    set(hObject,'Checked','on');
end

return

% --------------------------------------------------------------------
function plotMeanTseriesMenu_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function plotMeanTseriesCurrentCurrent_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
groupNum = viewGet(view,'currentGroup');
roiList = viewGet(view,'currentROI');
scanList = viewGet(view,'currentScan');
plotMeanTSeries(view, groupNum, roiList, scanList);

% --------------------------------------------------------------------
function plotMeanTseriesCurrentSelect_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
groupNum = viewGet(view,'currentGroup');
roiList = viewGet(view,'currentROI');
scanList = selectInList(view,'scans');
if ~isempty(scanList)
   plotMeanTSeries(view, groupNum, roiList, scanList);
end

% --------------------------------------------------------------------
function plotMeanTseriesCurrentAll_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
groupNum = viewGet(view,'currentGroup');
roiList = viewGet(view,'currentROI');
scanList = [1:viewGet(view,'nscans')];
plotMeanTSeries(view, groupNum, roiList, scanList);

% --------------------------------------------------------------------
function plotMeanTseriesAllCurrent_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
groupNum = viewGet(view,'currentGroup');
roiList = [1:viewGet(view,'numberofROIs')];
scanList = viewGet(view,'currentScan');
plotMeanTSeries(view, groupNum, roiList, scanList);

% --------------------------------------------------------------------
function plotMeanTseriesAllSelect_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
groupNum = viewGet(view,'currentGroup');
roiList = [1:viewGet(view,'numberofROIs')];
scanList = selectInList(view,'scans');
if ~isempty(scanList)
   plotMeanTSeries(view, groupNum, roiList, scanList);
%    plotMeanTSeries(view, groupNum, roiList, scanList,'subtractMean','No');
end

% --------------------------------------------------------------------
function plotMeanTseriesAllAll_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
groupNum = viewGet(view,'currentGroup');
roiList = [1:viewGet(view,'numberofROIs')];
scanList = [1:viewGet(view,'nscans')];
plotMeanTSeries(view, groupNum, roiList, scanList);

% --------------------------------------------------------------------
function plotMeanFourierAmpMenu_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function plotMeanFourierAmpCurrentCurrent_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
groupNum = viewGet(view,'currentGroup');
roiList = viewGet(view,'currentROI');
scanList = viewGet(view,'currentScan');
plotMeanFourierAmp(view, groupNum, roiList, scanList,  'detrend', 'Linear');

% --------------------------------------------------------------------
function plotMeanFourierAmpCurrentSelect_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
groupNum = viewGet(view,'currentGroup');
roiList = viewGet(view,'currentROI');
scanList = selectInList(view,'scans');
plotMeanFourierAmp(view, groupNum, roiList, scanList, 'detrend', 'Linear');

% --------------------------------------------------------------------
function plotMeanFourierAmpCurrentAll_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
groupNum = viewGet(view,'currentGroup');
roiList = viewGet(view,'currentROI');
scanList = [1:viewGet(view,'nscans')];
plotMeanFourierAmp(view, groupNum, roiList, scanList, 'detrend', 'Linear');

% --------------------------------------------------------------------
function plotMeanFourierAmpAllCurrent_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
groupNum = viewGet(view,'currentGroup');
roiList = [1:viewGet(view,'numberofROIs')];
scanList = viewGet(view,'currentScan');
plotMeanFourierAmp(view, groupNum, roiList, scanList, 'detrend', 'Linear');

% --------------------------------------------------------------------
function plotMeanFourierAmpAllSelect_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
groupNum = viewGet(view,'currentGroup');
roiList = [1:viewGet(view,'numberofROIs')];
scanList = selectInList(view,'scans');
plotMeanFourierAmp(view, groupNum, roiList, scanList, 'detrend', 'Linear');

% --------------------------------------------------------------------
function plotMeanFourierAmpAllAll_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
groupNum = viewGet(view,'currentGroup');
roiList = [1:viewGet(view,'numberofROIs')];
scanList = [1:viewGet(view,'nscans')];
plotMeanFourierAmp(view, groupNum, roiList, scanList, 'detrend', 'Linear');

% --------------------------------------------------------------------
function plotMotionCorrectionMatrices_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
d.groupNum = viewGet(view,'currentGroup');
d.scanNum = viewGet(view,'curScan');
d.expname = '';
d.ver = 4;
dispmotioncorrect(d);


% --------------------------------------------------------------------
function plotsDisplayEPIImages_Callback(hObject, eventdata, handles)

mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
mlrDisplayEPI(view);


% --------------------------------------------------------------------
function plotSpikeDetection_Callback(hObject, eventdata, handles)

mrGlobals;
viewNum = handles.viewNum;
view = MLR.views{viewNum};
mlrSpikeDetector(view,viewGet(view,'curScan'),viewGet(view,'curGroup'));


% --------------------------------------------------------------------
function flatViewerMenuItem_Callback(hObject, eventdata, handles)

mrGlobals;
viewNum = handles.viewNum;
v = MLR.views{viewNum};

% get the flat files info which
% is stored in the baseCoordMap
params = viewGet(v,'baseCoordMap');
flatPath = viewGet(v,'baseCoordMapPath');

% get what type of base this is (e.g. surface or flat)
baseType = viewGet(v,'baseType');
if baseType == 1
  % get the filename
  filename = sprintf('%s.off', stripext(params.flatFileName));
else
  % get the filename
  filename = sprintf('%s.off', stripext(params.outerSurfaceFileName));
end

% switch directories to the flatDir, asking
% the user to find it if it does not exist
thispwd = pwd;
if isdir(flatPath) && isfile(fullfile(flatPath,filename))
  cd(flatPath);
else
  if baseType ~= 1
    mrWarnDlg(sprintf('(mrLoadRetGUI) File %s does not exist, please find the anatomy folder for this surface/flat',fullfile(flatPath,filename)));
    pathStr = uigetdir(mrGetPref('volumeDirectory'),'Find anatomy folder from which this flat was created');
    if pathStr == 0,return,end
    cd(pathStr);
    viewSet(v,'baseCoordMapPath',pathStr);
  else
    % try to create the off for a flat
    flatParentSurf = fullfile(params.path,params.innerCoordsFileName);
    if isfile(flatParentSurf)
      disp('(mrLoadRetGUI) Creating missing flat off surface');
      disppercent(-inf,sprintf('(mrLoadRetGUI) Note this will create a quick flat surface good enough for rough visualization of location but is not exactly correct'));
      % load the parent surface 
      flatParentSurfOFF = loadSurfOFF(flatParentSurf);
      if ~isempty(flatParentSurfOFF)
	% find all vtcs that intersect
	flatVtcs = [];
	for corticalDepth = -0.5:0.05:1.5
	  flatVtcs = [flatVtcs ; reshape(params.innerCoords + (params.outerCoords-params.innerCoords)*corticalDepth,size(params.innerCoords,1)*size(params.innerCoords,2),3)];
	end
	% match coordinate to find vertices that this was made from
	%vtcs = assignToNearest(flatVtcs,flatParentSurfOFF.vtcs);
	% round and linearize the vertices
	flatVtcs = round(flatVtcs);
	flatLinear = mrSub2ind(params.dims,flatVtcs(:,1), flatVtcs(:,2),flatVtcs(:,3));
	parentVtcs = round(flatParentSurfOFF.vtcs);
	flatParentLinear = mrSub2ind(params.dims,parentVtcs(:,1),parentVtcs(:,2),parentVtcs(:,3));
	% find ones that match - note that this is an approximation
	% by finding vertices that are wihtin a rounded mm to
	% each other - also above we add across cortical depths.
	% Could be more precise - but that might be slower
	% and for this purposes (visualization of roughly
	% where the surface is, seems good enough)
	vtcs = find(ismember(flatParentLinear,flatLinear));
	% find tris that these vtcs are associated with
	matchingTris = flatParentSurfOFF.tris;
	[triNum edgeNum] = find(ismember(matchingTris,vtcs));
	tris = flatParentSurfOFF.tris(triNum,:);
	% now reget the vertices that are in all the triangles
	vtcs = unique(tris(:));
	% and renumber the vtcs in the tris 
	[dummy tris(:)] = ismember(tris(:),vtcs);
	% ok, put it all togehter
	flatSurf.filename = '';
	flatSurf.nParent = [flatParentSurfOFF.Nvtcs flatParentSurfOFF.Ntris 1]; 
	flatSurf.nPatch = [length(vtcs) size(tris,1) 1];
	flatSurf.vtcs = flatParentSurfOFF.vtcs(vtcs,:);
	flatSurf.tris = tris;
	flatSurf.parentSurfaceName = flatParentSurf;
	flatSurf.Nvtcs = length(vtcs);
	flatSurf.Ntris = size(tris,1);
	flatSurf.edges = 0;
	flatSurf.patch2parent(:,1) = (1:length(vtcs))';
	flatSurf.patch2parent(:,2) = vtcs';
	flatSurf.path = params.path;
	% put it into the params field
	params.flatFileName = flatSurf;
	disppercent(inf);
      end
    end
  end
end

if baseType == 1
  % now bring up the flat viewer
  mrFlatViewer(params.flatFileName,params.outerCoordsFileName,params.innerCoordsFileName,params.curvFileName,params.anatFileName,viewNum);
else
  mrSurfViewer(params.outerSurfaceFileName,params.outerCoordsFileName,params.innerSurfaceFileName,params.innerCoordsFileName,params.curvFileName,params.anatFileName);
end

cd(thispwd);

% --------------------------------------------------------------------
function calcDistMenu_Callback(hObject, eventdata, handles)
% hObject    handle to flatViewerMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function calcDistMenuItemSeg_Callback(hObject, eventdata, handles)
mrGlobals;
viewNum = handles.viewNum;
v = MLR.views{viewNum};
% [view userCancel] = newROI(view);
% if userCancel,return,end
if ~exist('dijkstra') == 3
  mrWarnDlg('(mrLoadRetGUI) You need to install the mrGray tools to use the function calcDist');
else
  [dijkstraDistance,euclidianDistance] = calcDist(v, 'segments'); % ,'polygon',1);
  if length(dijkstraDistance)==1
      disp(sprintf('Dijkstra Distance betwen pts A and B on surface: %2.4f', dijkstraDistance));
      disp(sprintf('Euclidian 3D distance betwen pts A and B: %2.4f', euclidianDistance));
  elseif length(dijkstraDistance)>1
    disp('Dijkstra Distance between successive pairs of points on surface:');
    disp(dijkstraDistance');
    disp('Euclidian 3D distance betwen successive pairs of points');
    disp(euclidianDistance');
  end
  refreshMLRDisplay(viewNum);
end


% --------------------------------------------------------------------
function editStimfileMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to editStimfileMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
mrGlobals;
viewNum = handles.viewNum;
v = MLR.views{viewNum};

% check for mgl function editStimfile
if exist('editStimfile') ~= 2
  mrWarnDlg(sprintf('(mrLoadRetGUI) To use this functionality, you need to have mgl installed. See http://justingardner.net/mgl'));
else
  editStimfile(v);
end


% --------------------------------------------------------------------
function scanViewInMlrVolMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to scanViewInMlrVolMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
mrGlobals;
viewNum = handles.viewNum;
v = MLR.views{viewNum};
mlrVol(v);
