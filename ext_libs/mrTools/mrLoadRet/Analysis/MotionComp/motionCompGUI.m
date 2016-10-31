function varargout = motionCompGUI(varargin)
% params = motionCompGUI('groupName',groupName,'params',params)
% params = motionCompGUI('groupName',groupName);
% params = motionCompGUI('params',params);
%
% GUI for averageTSeries.
% This function was created along with motionCompGUI.fig using GUIDE.
%
% params: parameter structure as in averageTSeries
%
% Examples:
%
% params = motionCompGUI('groupName','Raw');
% params = motionCompGUI('params',params);
% params = motionCompGUI('groupName','Raw','params',params);
%
% djh, 7/2004

% Last Modified by GUIDE v2.5 31-Jul-2007 13:26:41

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
  'gui_Singleton',  gui_Singleton, ...
  'gui_OpeningFcn', @motionCompGUI_OpeningFcn, ...
  'gui_OutputFcn',  @motionCompGUI_OutputFcn, ...
  'gui_LayoutFcn',  [] , ...
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


% --- Executes just before motionCompGUI is made visible.
function motionCompGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to motionCompGUI (see VARARGIN)

% Parse varargin
for index = 1:2:length(varargin)
  field = varargin{index};
  val = varargin{index+1};
  switch field
    case 'groupName'
      groupName = val;
    case 'params'
      params = val;
    otherwise
      mrWarnDlg('Invalid initialization argument')
  end
end

% Error if neither params nor groupName specified
if ieNotDefined('params') & ieNotDefined('groupName')
  mrErrorDlg('Must initialize with either params or groupName.');
end

% If groupName not passed then set it according to params.groupName (which
% we now know must exist when groupName does not).
if ieNotDefined('groupName')
  if ~isfield(params,'groupName')
    mrErrorDlg('Must initialize with params.groupName.');
  end
  groupName = params.groupName;
end

% Group names, group number, and nScans
groupNames = viewGet([],'groupNames');
if ~any(strcmp('MotionComp',groupNames))
  groupNames{length(groupNames)+1} = 'MotionComp';
end
groupNum = viewGet([],'groupNum',groupName);
if isempty(groupNum)
  mrErrorDlg('group ',groupName,' not found.');
end
nScans = viewGet([],'nscans',groupNum);

% Reconcile/initialize params with current status of group and ensure that
% it has the required fields.
if ieNotDefined('params')
  params = motionCompReconcileParams(groupName);
else
  params = motionCompReconcileParams(groupName,params);
end

% Initialize handles. Store the parameters in the GUI until user
% clicks the ok or cancel buttons. Then return the params from
% motionCompGUI_OutputFcn.
handles.groupName = params.groupName;
handles.crop = params.crop;

% Initialize include (targetScans)
handles.include = zeros(1,nScans);
handles.include(params.targetScans) = 1;

% Initialize tseriesfiles
handles.tseriesfiles = cell(1,nScans);
for scan = 1:nScans
  handles.tseriesfiles{scan} = viewGet([],'tseriesFile',scan,groupNum);
end

% Initialize descriptions
handles.descriptions = cell(1,nScans);
for scan = 1:nScans
  m = find(scan == params.targetScans);
  if m
    handles.descriptions{scan} = params.descriptions{m};
  else
    handles.descriptions{scan} = ['Motion compensation of ',groupName,' scan ',num2str(scan)];
  end
end

% Initialize group popup
groupNames = {'New',groupNames{:}};
motionCompGroupName = params.motionCompGroupName;
motionCompGroupNum = find(strcmp(motionCompGroupName,groupNames));
if isempty(motionCompGroupNum)
  groupNames{length(groupNames)+1} = motionCompGroupName;
  motionCompGroupNum = length(groupNames);
end
set(handles.motionCompGroupPopup,'String',groupNames);
set(handles.motionCompGroupPopup,'Value',motionCompGroupNum);

% Initialize base scan popup
baseScan = params.baseScan;
baseScanStrings = cell(nScans,1);
for scan = 1:nScans
  baseScanStrings{scan} = num2str(scan);
end
set(handles.baseScanPopup,'String',baseScanStrings);
set(handles.baseScanPopup,'Value',baseScan);

% Initialize base frame popup
baseFrameStrings = get(handles.baseFramePopup,'String');
baseFrameVal = find(strcmp(params.baseFrame,baseFrameStrings));
set(handles.baseFramePopup,'Value',baseFrameVal);

% Initialize interp method popup
interpMethodStrings = get(handles.interpMethodPopup,'String');
interpMethodVal = find(strcmp(params.interpMethod,interpMethodStrings));
set(handles.interpMethodPopup,'Value',interpMethodVal);

% Initialize checkboxes
set(handles.sliceTimeCheckbox,'Value',params.sliceTimeCorrection);
set(handles.robustCheckbox,'Value',params.robust);
set(handles.intensityContrastCheckbox,'Value',params.correctIntensityContrast);

% Initialize niters
set(handles.nitersText,'String',num2str(params.niters));

% Initialize current scan
setScan(handles,1);

% Initialize cancel flag
handles.cancel = 0;

% Update handles structure to include params and cancel
guidata(hObject, handles);

% UIWAIT makes motionCompGUI wait for user response (see UIRESUME)
% uiresume is called when user pushes ok or cancel buttons.
% This triggers the output function to return the params.
uiwait(handles.motionCompFigure);

% --- Outputs from this function are returned to the command line.
function varargout = motionCompGUI_OutputFcn(hObject, eventdata, handles)
% Return params or [] (if cancelled)
if handles.cancel
  varargout{1} = [];
else
  % Set params from the GUI
  params.groupName = handles.groupName;
  motionCompGroupNames = get(handles.motionCompGroupPopup,'String');
  motionCompGroupValue = get(handles.motionCompGroupPopup,'Value');
  params.motionCompGroupName = motionCompGroupNames{motionCompGroupValue};
  interpMethodStrings = get(handles.interpMethodPopup,'String');
  interpMethodValue = get(handles.interpMethodPopup,'Value');
  params.interpMethod = interpMethodStrings{interpMethodValue};
  params.baseScan = get(handles.baseScanPopup,'Value');
  baseFrameStrings = get(handles.baseFramePopup,'String');
  baseFrameValue = get(handles.baseFramePopup,'Value');
  params.baseFrame = baseFrameStrings{baseFrameValue};
  sliceTimeStrings = get(handles.sliceTimePopup,'String');
  sliceTimeValue = get(handles.sliceTimePopup,'Value');
  params.sliceTimeString = sliceTimeStrings{sliceTimeValue};
  params.sliceTimeCorrection = get(handles.sliceTimeCheckbox,'Value');
  params.robust = get(handles.robustCheckbox,'Value');
  params.correctIntensityContrast = get(handles.intensityContrastCheckbox,'Value');
  params.crop = handles.crop;
  params.niters = str2double(get(handles.nitersText,'String'));
  targetScans = find(handles.include);
  params.targetScans = targetScans;
  params.tseriesfiles = handles.tseriesfiles(targetScans);
  params.descriptions = handles.descriptions(targetScans);
  % Reconcile params and return
  params = motionCompReconcileParams(params.groupName,params);
  varargout{1} = params;
end
% Finally, close the GUI
close(handles.motionCompFigure);

% --- forwardButton.
function forwardButton_Callback(hObject, eventdata, handles)
curScan = str2double(get(handles.scanText,'String'));
setScan(handles,curScan+1);

% --- backwardButton.
function backwardButton_Callback(hObject, eventdata, handles)
curScan = str2double(get(handles.scanText,'String'));
setScan(handles,curScan-1);

% --- scanText
function scanText_CreateFcn(hObject, eventdata, handles)
if ispc
  set(hObject,'BackgroundColor','white');
else
  set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function scanText_Callback(hObject, eventdata, handles)
curScan = str2double(get(handles.scanText,'String'));
setScan(handles,curScan);

function setScan(handles,curScan)
% Update GUI to reflect curScan
% scan must be an integer between 1 and nScans
curScan = round(curScan);
nScans = length(handles.include);
curScan = min(max(1,curScan),nScans);
% scan text
set(handles.scanText,'String',num2str(curScan));
% description
groupNum = viewGet([],'groupNum',handles.groupName);
description = viewGet([],'description',curScan,groupNum);
set(handles.descriptionText,'String',description);
% include
set(handles.includeCheckbox,'Value',handles.include(curScan));

% --- includeCheckbox.
function includeCheckbox_Callback(hObject, eventdata, handles)
curScan = str2double(get(handles.scanText,'String'));
value = get(hObject,'Value');
handles.include(curScan) = value;
guidata(hObject, handles);

% --- sliceTimeCheckbox.
function sliceTimeCheckbox_Callback(hObject, eventdata, handles)
% value = get(hObject,'Value');
% handles.sliceTimeCorrection = value;
% guidata(hObject, handles);

% --- robustCheckbox.
function robustCheckbox_Callback(hObject, eventdata, handles)
% value = get(hObject,'Value');
% handles.robust = value;
% guidata(hObject, handles);

% --- intensityContrastCheckbox.
function intensityContrastCheckbox_Callback(hObject, eventdata, handles)
% value = get(hObject,'Value');
% handles.correctIntensityContrast = value;
% guidata(hObject, handles);

% --- nitersText
function nitersText_Callback(hObject, eventdata, handles)
% curScan = str2double(get(handles.scanText,'String'));
% value = round(str2double(get(hObject,'String')));
% set(hObject,'String',num2str(value));
% handles.niters(curScan) = value;
% guidata(hObject, handles);

function nitersText_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
  set(hObject,'BackgroundColor','white');
end

% --- motionCompGroupPopup.
function motionCompGroupPopup_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
  set(hObject,'BackgroundColor','white');
end

function motionCompGroupPopup_Callback(hObject, eventdata, handles)
strings = get(hObject,'String');
value = get(hObject,'Value');
if (value == 1)
  options.Resize = 'off';
  options.WindowStyle = 'modal';
  options.Interpreter = 'none';
  newGroupName = inputdlg('Enter name for new group: ','New group',1,{''},options);
  if ~isempty(newGroupName)
    strings = {strings{:},newGroupName{1}};
    set(hObject,'String',strings);
    set(hObject,'Value',length(strings));
  end
end
guidata(hObject, handles);

% --- baseScanPopup.
function baseScanPopup_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
  set(hObject,'BackgroundColor','white');
end

function baseScanPopup_Callback(hObject, eventdata, handles)

% --- interpMethodPopup.
function interpMethodPopup_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
  set(hObject,'BackgroundColor','white');
end

function interpMethodPopup_Callback(hObject, eventdata, handles)

% --- baseFramePopup.
function baseFramePopup_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
  set(hObject,'BackgroundColor','white');
end

function baseFramePopup_Callback(hObject, eventdata, handles)

% --- sliceTimePopup.
function sliceTimePopup_Callback(hObject, eventdata, handles)

function sliceTimePopup_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- cropButton.
function cropButton_Callback(hObject, eventdata, handles)
view = newView;
baseScan = get(handles.baseScanPopup,'Value');
% Load first frame
volume = loadTSeries(view,baseScan,'all',1);
handles.crop = selectCropRegion(volume);
% delete temporary view
deleteView(view);
guidata(hObject, handles);

% --- okButton.
function okButton_Callback(hObject, eventdata, handles)
% uiresume triggers motionCompGUI_OutputFcn
uiresume;

% --- cancelButton.
function cancelButton_Callback(hObject, eventdata, handles)
handles.cancel = 1;
guidata(hObject, handles);
% uiresume triggers motionCompGUI_OutputFcn
uiresume;
return




