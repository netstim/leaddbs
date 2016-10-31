function varargout = corAnalGUI(varargin)
% params = corAnalGUI('groupName',groupName,'params',params)
% params = corAnalGUI('groupName',groupName);
% params = corAnalGUI('params',params);
% 
% Creates a new corAnal GUI.
% This function was created along with mrLoadRetGui.fig using GUIDE.
%
% params: parameter structure as in averageTSeries
%
% Examples:
%
% params = corAnalGUI('groupName','Raw');
% params = corAnalGUI('params',params);
% params = corAnalGUI('groupName','Raw','params',params);
%
% djh, 7/2004

% Last Modified by GUIDE v2.5 15-Nov-2010 16:59:48

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @corAnalGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @corAnalGUI_OutputFcn, ...
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


% --- Executes just before corAnalGUI is made visible.
function corAnalGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to corAnalGUI (see VARARGIN)

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

% Reconcile/initialize params with current status of group and ensure that
% params has the required fields.
if ieNotDefined('params')
    params = corAnalReconcileParams(groupName);
else
    params = corAnalReconcileParams(groupName,params);
end

% Group names, group number, and nScans
groupNames = viewGet([],'groupNames');
groupNum = viewGet([],'groupNum',groupName);
if isempty(groupNum)
    mrErrorDlg('group ',groupName,' not found.');
end
nScans = viewGet([],'nscans',groupNum);

% Initialize handles. Store the parameters in the GUI until user
% clicks the ok or cancel buttons. Then return the params from
% corAnalGUI_OutputFcn.
handles.groupName = params.groupName;
handles.recompute = params.recompute;
handles.ncycles = params.ncycles;
handles.detrend = params.detrend;
handles.spatialnorm = params.spatialnorm;
handles.trigonometricFunction = params.trigonometricFunction;
handles.tseriesfile = params.tseriesfile;
setScan(handles,1);

% Initialize cancel flag
handles.cancel = 0;

% Update handles structure to include params and cancel
guidata(hObject, handles);

% UIWAIT makes corAnalGUI wait for user response (see UIRESUME)
% uiresume is called when user pushes ok or cancel buttons.
% This triggers the output function to return the params.
uiwait(handles.corAnalFigure);


% --- Outputs from this corAnalGUI 
function varargout = corAnalGUI_OutputFcn(hObject, eventdata, handles)
% Return params or [] (if cancelled)
if handles.cancel
    varargout{1} = [];
else
    % Set params from the GUI
    params.groupName = handles.groupName;
    params.recompute = handles.recompute;
    params.ncycles = handles.ncycles;
    params.detrend = handles.detrend;
    params.spatialnorm = handles.spatialnorm;
    params.trigonometricFunction = handles.trigonometricFunction;
    params.tseriesfile = handles.tseriesfile;
    % Reconcile params and return
    params = corAnalReconcileParams(params.groupName,params);
    varargout{1} = params;
end
% Finally, close the GUI
close(handles.corAnalFigure);

% --- forwardButton
function forwardButton_Callback(hObject, eventdata, handles)
curScan = str2double(get(handles.scanText,'String'));
setScan(handles,curScan+1);

% --- backwardButton
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
nScans = length(handles.recompute);
curScan = min(max(1,curScan),nScans);
% scan text
set(handles.scanText,'String',num2str(curScan));  
% description
groupNum = viewGet([],'groupNum',handles.groupName);
description = viewGet([],'description',curScan,groupNum);
set(handles.descriptionText,'String',description);
% recompute
set(handles.recomputeCheckbox,'Value',handles.recompute(curScan));
% ncycles
set(handles.nCyclesText,'String',num2str(handles.ncycles(curScan)));
% detrend
detrendValue = find(strcmp(handles.detrend{curScan},...
    get(handles.detrendPopup,'String')));
set(handles.detrendPopup,'Value',detrendValue);
% spatialnormlabel
spatialnormValue = find(strcmp(handles.spatialnorm{curScan},...
    get(handles.spatialNormPopup,'String')));
set(handles.spatialNormPopup,'Value',spatialnormValue);
% trigonometricFunctionlabel
trigonometricFunctionValue = find(strcmp(handles.trigonometricFunction{curScan},...
    get(handles.trigonometricFunctionPopup,'String')));
set(handles.trigonometricFunctionPopup,'Value',trigonometricFunctionValue);

% --- copyButton.
function copyButton_Callback(hObject, eventdata, handles)
curScan = str2double(get(handles.scanText,'String'));
nScans = length(handles.recompute);
for scan = 1:nScans
    handles.recompute(scan) = handles.recompute(curScan);
    handles.ncycles(scan) = handles.ncycles(curScan);
    handles.detrend{scan} = handles.detrend{curScan};
    handles.spatialnorm{scan} = handles.spatialnorm{curScan};
    handles.trigonometricFunction{scan} = handles.trigonometricFunction{curScan};
end
guidata(hObject, handles);

% --- recomputeCheckbox.
function recomputeCheckbox_Callback(hObject, eventdata, handles)
% Hint: get(hObject,'Value') returns toggle state of recomputeCheckbox
curScan = str2double(get(handles.scanText,'String'));
value = get(hObject,'Value');
handles.recompute(curScan) = value;
guidata(hObject, handles);

% --- nCyclesText
function nCyclesText_Callback(hObject, eventdata, handles)
curScan = str2double(get(handles.scanText,'String'));
value = round(str2double(get(hObject,'String')));
set(hObject,'String',num2str(value)); 
handles.ncycles(curScan) = value;
guidata(hObject, handles);

function nCyclesText_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- detrendPopup
function detrendPopup_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function detrendPopup_Callback(hObject, eventdata, handles)
% Hints: contents = get(hObject,'String') returns detrendPopup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from detrendPopup
curScan = str2double(get(handles.scanText,'String'));
string = get(hObject,'String');
value = get(hObject,'Value');
handles.detrend{curScan} = string{value};
guidata(hObject, handles);

% --- spatialNormPopup
function spatialNormPopup_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function spatialNormPopup_Callback(hObject, eventdata, handles)
% Hints: contents = get(hObject,'String') returns spatialNormPopup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from spatialNormPopup
curScan = str2double(get(handles.scanText,'String'));
string = get(hObject,'String');
value = get(hObject,'Value');
handles.spatialnorm{curScan} = string{value};
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function trigonometricFunctionPopup_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in trigonometricFunctionPopup.
function trigonometricFunctionPopup_Callback(hObject, eventdata, handles)
curScan = str2double(get(handles.scanText,'String'));
string = get(hObject,'String');
value = get(hObject,'Value');
handles.trigonometricFunction{curScan} = string{value};
guidata(hObject, handles);

% --- okButton.
function okButton_Callback(hObject, eventdata, handles)
% uiresume triggers corAnalGUI_OutputFcn
uiresume;

% --- cancelButton.
function cancelButton_Callback(hObject, eventdata, handles)
handles.cancel = 1;
guidata(hObject, handles);
% uiresume triggers corAnalGUI_OutputFcn
uiresume;
return




