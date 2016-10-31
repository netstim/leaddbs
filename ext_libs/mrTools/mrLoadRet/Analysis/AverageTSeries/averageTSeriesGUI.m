function varargout = averageTSeriesGUI(varargin)
% params = averageTSeriesGUI('groupName',groupName,'params',params)
% params = averageTSeriesGUI('groupName',groupName);
% params = averageTSeriesGUI('params',params);
% 
% GUI for averageTSeries.
% This function was created along with averageTSeriesGUI.fig using GUIDE.
%
% params: parameter structure as in averageTSeries
%
% Examples:
%
% params = averageTSeriesGUI('groupName','Raw');
% params = averageTSeriesGUI('params',params);
% params = averageTSeriesGUI('groupName','Raw','params',params);
%
% djh, 7/2004

% Last Modified by GUIDE v2.5 26-Jul-2006 17:59:37

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @averageTSeriesGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @averageTSeriesGUI_OutputFcn, ...
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


% --- Executes just before averageTSeriesGUI is made visible.
function averageTSeriesGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to averageTSeriesGUI (see VARARGIN)

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
if ~any(strcmp('Averages',groupNames))
    groupNames{length(groupNames)+1} = 'Averages';
end
groupNum = viewGet([],'groupNum',groupName);
if isempty(groupNum)
    mrErrorDlg('group ',groupName,' not found.');
end
nScans = viewGet([],'nscans',groupNum);

% Reconcile/initialize params with current status of group and ensure that
% it has the required fields.
if ieNotDefined('params')
    params = averageTSeriesReconcileParams(groupName);
else
    params = averageTSeriesReconcileParams(groupName,params)
end

% Initialize handles. Store the parameters in the GUI until user
% clicks the ok or cancel buttons. Then return the params from
% averageTSeriesGUI_OutputFcn.
handles.groupName = params.groupName;
handles.description = params.description;
handles.fileName = [];

% Initialize tseriesfiles
handles.tseriesfiles = cell(1,nScans);
for scan = 1:nScans
    handles.tseriesfiles{scan} = viewGet([],'tseriesFile',scan,groupNum);
end

% Initialize include (scanList)
handles.include = zeros(1,nScans);
handles.include(params.scanList) = 1;

% Initalize shift and revers
handles.shift = zeros(1,nScans);
handles.shift(params.scanList) = params.shiftList;
handles.reverse = zeros(1,nScans);
handles.reverse(params.scanList) = params.reverseList;

% Initialize interp method popup
interpMethodStrings = get(handles.interpMethodPopup,'String');
interpMethodVal = find(strcmp(params.interpMethod,interpMethodStrings));
set(handles.interpMethodPopup,'Value',interpMethodVal);

% Initialize base scan popup
baseScan = params.baseScan;
baseScanStrings = cell(nScans,1);
for scan = 1:nScans
	baseScanStrings{scan} = num2str(scan);
end
set(handles.baseScanPopup,'String',baseScanStrings);
set(handles.baseScanPopup,'Value',baseScan);

% Initialize aveGroup popup
groupNames = {'New',groupNames{:}};
aveGroupName = params.aveGroupName;
aveGroupNum = find(strcmp(aveGroupName,groupNames));
if isempty(aveGroupNum)
    groupNames{length(groupNames)+1} = aveGroupName;
    aveGroupNum = length(groupNames);
end
set(handles.aveGroupPopup,'String',groupNames);
set(handles.aveGroupPopup,'Value',aveGroupNum);

% Initialize current scan
setScan(handles,1);

% Initialize cancel flag
handles.cancel = 0;

% Update handles structure to include params and cancel
guidata(hObject, handles);

% UIWAIT makes averageTSeriesGUI wait for user response (see UIRESUME)
% uiresume is called when user pushes ok or cancel buttons.
% This triggers the output function to return the params.
uiwait(handles.averageTSeriesFigure);

% --- Outputs from this function are returned to the command line.
function varargout = averageTSeriesGUI_OutputFcn(hObject, eventdata, handles)
% Return params or [] (if cancelled)
if handles.cancel
    varargout{1} = [];
else
    % Set params from the GUI
	scanList = find(handles.include);
	params.scanList = scanList;
    params.tseriesfiles = handles.tseriesfiles(scanList);
	params.shiftList = handles.shift(scanList);
	params.reverseList = handles.reverse(scanList);
	params.baseScan = get(handles.baseScanPopup,'Value');
	params.groupName = handles.groupName;
	aveGroupNames = get(handles.aveGroupPopup,'String');
	aveGroupValue = get(handles.aveGroupPopup,'Value');
	params.aveGroupName = aveGroupNames{aveGroupValue};
	interpMethodStrings = get(handles.interpMethodPopup,'String');
    interpMethodValue = get(handles.interpMethodPopup,'Value');
	params.interpMethod = interpMethodStrings{interpMethodValue};
	params.description = handles.description;
	params.fileName = handles.fileName;
    % Reconcile params and return
    params = averageTSeriesReconcileParams(params.groupName,params);
    varargout{1} = params;
end
% Finally, close the GUI
close(handles.averageTSeriesFigure);

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
% shift
set(handles.shiftText,'String',num2str(handles.shift(curScan)));
% reverse
set(handles.reverseCheckbox,'Value',handles.reverse(curScan));

% --- copyButton.
function copyButton_Callback(hObject, eventdata, handles)
curScan = str2double(get(handles.scanText,'String'));
nScans = length(handles.include);
for scan = 1:nScans
    handles.include(scan) = handles.include(curScan);
	handles.shift(scan) = handles.shift(curScan);
	handles.reverse(scan) = handles.reverse(curScan);
end
guidata(hObject, handles);

% --- includeCheckbox.
function includeCheckbox_Callback(hObject, eventdata, handles)
curScan = str2double(get(handles.scanText,'String'));
value = get(hObject,'Value');
handles.include(curScan) = value;
guidata(hObject, handles);

% --- reverseCheckbox.
function reverseCheckbox_Callback(hObject, eventdata, handles)
curScan = str2double(get(handles.scanText,'String'));
value = get(hObject,'Value');
handles.reverse(curScan) = value;
guidata(hObject, handles);

% --- shiftText
function shiftText_Callback(hObject, eventdata, handles)
curScan = str2double(get(handles.scanText,'String'));
value = round(str2double(get(hObject,'String')));
set(hObject,'String',num2str(value)); 
handles.shift(curScan) = value;
guidata(hObject, handles);

function shiftText_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- aveGroupPopup.
function aveGroupPopup_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function aveGroupPopup_Callback(hObject, eventdata, handles)
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

% --- interpMethodPopup.
function interpMethodPopup_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function interpMethodPopup_Callback(hObject, eventdata, handles)

% --- baseScanPopup.
function baseScanPopup_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function baseScanPopup_Callback(hObject, eventdata, handles)
guidata(hObject, handles);

% --- okButton.
function okButton_Callback(hObject, eventdata, handles)
% uiresume triggers averageTSeriesGUI_OutputFcn
uiresume;

% --- cancelButton.
function cancelButton_Callback(hObject, eventdata, handles)
handles.cancel = 1;
guidata(hObject, handles);
% uiresume triggers averageTSeriesGUI_OutputFcn
uiresume;
return

