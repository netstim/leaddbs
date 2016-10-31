function varargout = editGroupGUI(varargin)
% Initialize/modify scanParams
% editGroupGUI can be called by passing an existing group
%     group = editGroupGUI('group',group);
% or by passing a groupname in which case the structure is built from
% scratch
%     group = editGroupGUI('groupname',groupname);

% Last Modified by GUIDE v2.5 14-May-2005 14:00:26

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @editGroupGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @editGroupGUI_OutputFcn, ...
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


% --- Executes just before editGroupGUI is made visible.
function editGroupGUI_OpeningFcn(hObject, eventdata, handles, varargin)

if length(varargin) == 2
    field = varargin{1};
    val = varargin{2};
    switch field
        case 'group'
            group = val;
        case 'groupname'
            group.name = val;
            group.scanParams = getScanParams(val);
            group.auxParams = [];
        otherwise
            warn('editGroupGUI: invalid initialization argument');
    end
end

% Update handles structure
handles.group = group;
guidata(hObject, handles);

% Initialize gui
setScan(handles,1);

% Choose default command line output for editGroupGUI
handles.output = handles.group;

% UIWAIT makes editGroupGUI wait for user response (see UIRESUME)
uiwait(handles.figure);


% --- Outputs from this function are returned to the command line.
function varargout = editGroupGUI_OutputFcn(hObject, eventdata, handles) 
% Returns group
if exist('handles')
    varargout{1} = handles.output;
    delete(handles.figure);
else
    varargout{1} = [];
end

% --- getScanParams
function scanParams = getScanParams(groupname)
% Reads nifti headers, and loads that information into the
% scanParams struct array

mrGlobals
  
scanParams = [];
tseriesDir = fullfile(MLR.homeDir,groupname,'TSeries');
if ~exist(tseriesDir, 'dir')
  mrErrorDlg(['No TSeries directory found in ',tseriesDir]);
end

% Get file names for *.img, *.hdr, and *.nii files
[nfilesImg, fileListImg] = countNiftiFiles(tseriesDir,'.img');
[nfilesHdr, fileListHdr] = countNiftiFiles(tseriesDir,'.hdr');
[nfilesNii, fileListNii] = countNiftiFiles(tseriesDir,'.nii');

% Check that each *.img file has a corresponding *.hdr file
if (nfilesImg ~= nfilesHdr)
    mrErrorDlg('Header file numbers (*.hdr) are inconsistent with data file numbers (*.img)');
end
for f = 1:nfilesImg
    [path,fname1,ext] = fileparts(fileListImg{f});
    [path,fname2,ext] = fileparts(fileListHdr{f});
    if ~strcmp(fname1,fname2)
        mrErrorDlg('Header file numbers (*.hdr) are inconsistent with data file numbers (*.img)');
    end
end

% Concatenate file lists
nfiles = nfilesImg + nfilesNii;
fileList = {fileListImg{:},fileListNii{:}};
headerFileList = {fileListHdr{:},fileListNii{:}};

% fileList{1}
% headerFileList{1}

% Check that there are tSeries files and that they have matching header files.
if ~nfiles
    mrErrorDlg(['No Tseries files found in ', tseriesDir]);
end

for iScan=1:nfiles
  name = fullfile(tseriesDir, headerFileList{iScan});
  hdr = mlrImageReadNiftiHeader(name);
  scanParams(iScan).dataSize = hdr.dim([2,3,4])';
  scanParams(iScan).description = 'description';
  scanParams(iScan).fileName = fileList{iScan};
  scanParams(iScan).fileType = 'Nifti';
  % scanParams(iScan).framePeriod = hdr.pixdim(5)/1000;
  disp('(editGroupGUI) checking time units in nifti file')
  niftiSpaceUnit = rem(hdr.xyzt_units, 8); 
  niftiTimeUnit = rem(hdr.xyzt_units-niftiSpaceUnit, 64);
  if niftiTimeUnit == 8 % seconds
    scanParams(iScan).framePeriod = hdr.pixdim(5)./1;
  elseif niftiTimeUnit == 16 % milliseconds
    scanParams(iScan).framePeriod = hdr.pixdim(5)./1000;
  elseif niftiTimeUnit == 32 % microseconds
    scanParams(iScan).framePeriod = hdr.pixdim(5)./10e6;
  end
  if strcmp(lower(mrGetPref('verbose')),'yes')
    % 8 -> 10^0, 16 -> 10^3, 32-> 10^6
    disp(sprintf('(viewSet) Timing. Pixdim(5) units: %d. Scaling by 10e%d',niftiTimeUnit, 3*(log2(niftiTimeUnit)-3)));
  end
  scanParams(iScan).junkFrames = 0;
  scanParams(iScan).nFrames = hdr.dim(5);
  scanParams(iScan).niftiHdr = hdr;
  scanParams(iScan).originalFileName{1} = fileList{iScan};
  scanParams(iScan).originalGroupName{1} = groupname;
  scanParams(iScan).totalFrames = hdr.dim(5);
  scanParams(iScan).voxelSize = hdr.pixdim([2,3,4])';
end  


% --- forwardButton
function forwardButton_Callback(hObject, eventdata, handles)
curScan = str2double(get(handles.scanText,'String'));
setScan(handles,curScan+1);


% --- backwardButton
function backwardButton_Callback(hObject, eventdata, handles)
curScan = str2double(get(handles.scanText,'String'));
setScan(handles,curScan-1);


% --- scanText
function scanText_Callback(hObject, eventdata, handles)
curScan = str2double(get(handles.scanText,'String'));
setScan(handles,curScan);

function scanText_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- setScan
function setScan(handles,curScan)
% Update GUI to reflect curScan
% curScan must be an integer between 1 and nScans
curScan = round(curScan);
curScan = min(max(1,curScan),length(handles.group.scanParams));
set(handles.scanText,'string',num2str(curScan));
scanParams = handles.group.scanParams(curScan);
% file name
set(handles.fileName,'string',scanParams.fileName);
% description
set(handles.descriptionText,'string',scanParams.description);
% frames
set(handles.totalFrames,'string',['Total frames: ',num2str(scanParams.totalFrames)]);
set(handles.nFramesText,'string',num2str(scanParams.nFrames));
set(handles.junkFramesText,'string',num2str(scanParams.junkFrames));


% --- copyButton
function copyButton_Callback(hObject, eventdata, handles)
curScan = str2double(get(handles.scanText,'String'));
nScans = length(handles.group.scanParams);
for scan = 1:nScans
    handles.group.scanParams(scan).description = handles.group.scanParams(curScan).description;
    totalFrames = handles.group.scanParams(scan).totalFrames;
    junkFrames = handles.group.scanParams(curScan).junkFrames;
    junkFrames = min([junkFrames,totalFrames]);
    handles.group.scanParams(scan).junkFrames = junkFrames;
    nFrames = totalFrames - junkFrames;
    handles.group.scanParams(scan).nFrames = nFrames;
end
guidata(hObject, handles);


% --- description
function descriptionText_Callback(hObject, eventdata, handles)
% Update group.scanParams.foo
curScan = str2double(get(handles.scanText,'string'));
handles.group.scanParams(curScan).description = get(hObject,'string');
guidata(hObject, handles);

function descriptionText_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- junkFrames
function junkFramesText_Callback(hObject, eventdata, handles)
junkFrames = round(str2double(get(hObject,'string')));
curScan = str2double(get(handles.scanText,'String'));
if isfinite(junkFrames)
    totalFrames = handles.group.scanParams(curScan).totalFrames;
    junkFrames = min([junkFrames,totalFrames]);
    handles.group.scanParams(curScan).junkFrames = junkFrames;
    handles.group.scanParams(curScan).nFrames = totalFrames - junkFrames;
    guidata(hObject, handles);
end
setScan(handles,curScan);

function junkFramesText_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- nFrames
function nFramesText_Callback(hObject, eventdata, handles)
nFrames = round(str2double(get(hObject,'string')));
curScan = str2double(get(handles.scanText,'String'));
if isfinite(nFrames)
    totalFrames = handles.group.scanParams(curScan).totalFrames;
    nFrames = min([nFrames,totalFrames]);
    handles.group.scanParams(curScan).nFrames = nFrames;
%    handles.group.scanParams(curScan).junkFrames = totalFrames - nFrames;
    guidata(hObject, handles);
end
setScan(handles,curScan);

function nFramesText_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- OK
function okButton_Callback(hObject, eventdata, handles)
% Check that required fields are filled
if ~isgroup(handles.group)
    mrErrorDlg('Invalid group');
end
handles.output = handles.group;
guidata(hObject, handles);
uiresume;


% --- Cancel
function cancelButton_Callback(hObject, eventdata, handles)
handles.output = [];
guidata(hObject, handles);
uiresume;




