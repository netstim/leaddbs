function varargout = mrAlignGUI(varargin)
% mrAlignGUI M-file for mrAlignGUI.fig
% See also: GUIDE, GUIDATA, GUIHANDLES
%        $Id$

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @mrAlignGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @mrAlignGUI_OutputFcn, ...
                   'gui_LayoutFcn',  [], ...
                   'gui_Callback',   []);
if nargin && isstr(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before mrAlignGUI is made visible.
function mrAlignGUI_OpeningFcn(hObject, eventdata, handles, varargin)
global ALIGN

% Initialize ALIGN global variable
ALIGN.volumePath = [];
ALIGN.inplanePath = [];
ALIGN.xform = eye(4);
ALIGN.guiXform = eye(4);
ALIGN.volSize = [64 64 64];
ALIGN.inplaneSize = [64 64 64];

% set the location of the figure
figloc = mrGetFigLoc('mrAlignGUI');
if ~isempty(figloc)
    set(handles.figure1,'Position',figloc);
end

% Initialize GUI
mrAlignGUI('sagittalRadioButton_Callback',hObject, eventdata, handles);
set(handles.transposeButton,'Value',0);
set(handles.flipButton,'Value',0);
set(handles.overlayButton,'Value',1);
set(handles.transparencySlider,'Value',1);
set(handles.setTalXform,'Enable','off');
set(handles.exportTal2Session,'Enable','off');
setAlignGUI(handles,'rot',[0 0 0]);
setAlignGUI(handles,'trans',[0 0 0]);
refreshAlignDisplay(handles);

% Choose default command line output for mrAlignGUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% --- Executes during object deletion, before destroying properties.
function axes_DeleteFcn(hObject, eventdata, handles)
clear global ALIGN

% --- Outputs from this function are returned to the command line.
function varargout = mrAlignGUI_OutputFcn(hObject, eventdata, handles)
% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on mouse press over axes background.
function axes_ButtonDownFcn(hObject, eventdata, handles)

% --- Executes on button press in sagittalRadioButton.
function sagittalRadioButton_Callback(hObject, eventdata, handles)
global ALIGN
set(handles.sagittalRadioButton,'Value',1);
set(handles.coronalRadioButton,'Value',0);
set(handles.axialRadioButton,'Value',0);
[m,index] = max(ALIGN.volumePermutation * [1 0 0]');
ALIGN.sliceOrientation = index;
setAlignGUI(handles,'nSlices',ALIGN.volSize(ALIGN.sliceOrientation));
set(handles.sliceSlider,'value',ALIGN.coords(ALIGN.sliceOrientation));
refreshAlignDisplay(handles);

% --- Executes on button press in coronalRadioButton.
function coronalRadioButton_Callback(hObject, eventdata, handles)
global ALIGN
set(handles.sagittalRadioButton,'Value',0);
set(handles.coronalRadioButton,'Value',1);
set(handles.axialRadioButton,'Value',0);
[m,index] = max(ALIGN.volumePermutation * [0 1 0]');
ALIGN.sliceOrientation = index;
setAlignGUI(handles,'nSlices',ALIGN.volSize(ALIGN.sliceOrientation));
set(handles.sliceSlider,'value',ALIGN.coords(ALIGN.sliceOrientation));
refreshAlignDisplay(handles);

% --- Executes on button press in axialRadioButton.
function axialRadioButton_Callback(hObject, eventdata, handles)
global ALIGN
set(handles.sagittalRadioButton,'Value',0);
set(handles.coronalRadioButton,'Value',0);
set(handles.axialRadioButton,'Value',1);
[m,index] = max(ALIGN.volumePermutation * [0 0 1]');
ALIGN.sliceOrientation = index;
setAlignGUI(handles,'nSlices',ALIGN.volSize(ALIGN.sliceOrientation));
set(handles.sliceSlider,'value',ALIGN.coords(ALIGN.sliceOrientation));
refreshAlignDisplay(handles);

% --- Executes on button press in overlayButton.
function overlayButton_Callback(hObject, eventdata, handles)
refreshAlignDisplay(handles);

% --- Executes on button press in transposeButton.
function transposeButton_Callback(hObject, eventdata, handles)
refreshAlignDisplay(handles);

% --- Executes on button press in flipButton.
function flipButton_Callback(hObject, eventdata, handles)
refreshAlignDisplay(handles);

% --- Executes during object creation, after setting all properties.
function sliceSlider_CreateFcn(hObject, eventdata, handles)
usewhitebg = 1;
if usewhitebg
    set(hObject,'BackgroundColor',[.9 .9 .9]);
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

% --- Executes on slider movement.
function sliceSlider_Callback(hObject, eventdata, handles)
global ALIGN
slice = (round(get(hObject,'Value')));
ALIGN.coords(ALIGN.sliceOrientation) = slice;
refreshAlignDisplay(handles);

% --- Executes during object creation, after setting all properties.
function transparencySlider_CreateFcn(hObject, eventdata, handles)
usewhitebg = 1;
if usewhitebg
    set(hObject,'BackgroundColor',[.9 .9 .9]);
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

% --- Executes on slider movement.
function transparencySlider_Callback(hObject, eventdata, handles)
global ALIGN
refreshAlignDisplay(handles);

% --- Executes during object creation, after setting all properties.
function transX_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function transX_Callback(hObject, eventdata, handles)
global ALIGN
ALIGN.guiXform = getGuiXform(handles);
refreshAlignDisplay(handles);

% --- Executes during object creation, after setting all properties.
function transY_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function transY_Callback(hObject, eventdata, handles)
global ALIGN
ALIGN.guiXform = getGuiXform(handles);
refreshAlignDisplay(handles);

% --- Executes during object creation, after setting all properties.
function transZ_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function transZ_Callback(hObject, eventdata, handles)
global ALIGN
ALIGN.guiXform = getGuiXform(handles);
refreshAlignDisplay(handles);

% --- Executes during object creation, after setting all properties.
function rotX_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function rotX_Callback(hObject, eventdata, handles)
global ALIGN
ALIGN.guiXform = getGuiXform(handles);
refreshAlignDisplay(handles);

% --- Executes during object creation, after setting all properties.
function rotY_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function rotY_Callback(hObject, eventdata, handles)
global ALIGN
ALIGN.guiXform = getGuiXform(handles);
refreshAlignDisplay(handles);

% --- Executes during object creation, after setting all properties.
function rotZ_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function rotZ_Callback(hObject, eventdata, handles)
global ALIGN
ALIGN.guiXform = getGuiXform(handles);
refreshAlignDisplay(handles);


% --------------------------------------------------------------------
function fileMenu_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function loadVolMenuItem_Callback(hObject, eventdata, handles)

% Prompt user to choose volume.
initPath = mrGetPref('volumeDirectory');
if isempty(initPath)
    initPath = pwd;
end
pathStr = mlrGetPathStrDialog(initPath,'Choose vAnatomy file',{'*.img;*.nii','NIFTI Files'});
if isempty(pathStr),return,end

mrAlignLoadVol(pathStr,hObject,eventdata,handles);

function mrAlignLoadVol(pathStr,hObject,eventdata,handles)
% function load volume anatomy (taken out of loadVol) since
% it is reused by loadSourceAsDestination

global ALIGN
ALIGN.volumePath = pathStr;

% Load volume and header
h = mrMsgBox('Loading volume. Please wait');
[vData,hdr] = mlrImageReadNifti(ALIGN.volumePath);
mrCloseDlg(h);
if isempty(vData),return,end
volumeDimension = length(size(vData));

% Load associated .mat file for the volume (used in mrLR for bases,
% and needed here for keeping track of Talairach information). If
% no such .mat file exists, create it. (If user doesn't use mrLR,
% they can just ignore this extra file
matFilename = sprintf('%s.mat',stripext(pathStr));
if ~(exist(matFilename,'file')==2) % if base doesn't already have an associated mat file
  ALIGN.volBase.name = pathStr;
  ALIGN.volBase.data = vData;
  ALIGN.volBase.hdr = hdr;
  ALIGN.volBase.permutationMatrix = getPermutationMatrix(hdr);
  % check if it's a canonical volume, and if so, use qform as vol2mag
  roundErr = 100000;
  if ~isequal(hdr.sform_code,0) && isequal(round(hdr.qform44*roundErr),round(hdr.sform44*roundErr))
    mrWarnDlg(['Attention mrLoadRet users: This destination volume' ...
              ' seems to be a canonical base volume (The qform is' ...
	      ' equal to the sform). We are treating it ' ...
              ' as the main volume to which everything will be aligned.' ...
              ' In particular, we are setting its qform to be the vol2mag.' ...
              ' If this is correct, then you should save this by using Set' ...
	      ' Base Coordinate Frame in the File menu (and then you will' ...
	      ' not see this message again). If it is incorrect, you may' ...
	      ' want to reset the vol2mag.']');
    ALIGN.volBase.vol2mag = hdr.qform44;
  end % if not, will be automatically set to [] by isbase()
  % check if it's a talairach base, and if so, use sform as vol2tal
  if hdr.sform_code == 3 % if there's no .mat file but it's a Tal base
    mrWarnDlg(['Attention mrLoadRet users: This destination volume'...
              ' seems to be a canonical base volume with Talairach'...
              ' coordinates defined. Using the transformation of this volume' ...
              ' to Talairach space (e.g. its s-form) as the base.vol2tal,'...
              ' and using its transformation to magnet space (e.g. its'...
              ' qform) as base.vol2mag. If this is not correct, then you' ...
              ' need to fix the base structure in mrLoadRet.']);
    ALIGN.volBase.vol2tal = hdr.sform44;
    ALIGN.volBase.vol2mag = hdr.qform44;
    set(handles.exportTal2Session,'Enable','on');% allow user to export the vol2tal without redefining
  end % if not, will autmatically set vol2tal = []  
else % if there already is a base file
  load(matFilename); % then load it
  ALIGN.volBase = base; clear base;
  % if there is already a vol2tal, allow user to export it without redefining it
  if ~isempty(ALIGN.volBase.vol2tal), set(handles.exportTal2Session,'Enable','on'); end 
end
[tf ALIGN.volBase] = isbase(ALIGN.volBase); % make sure it has all the right fields;


% check rank of qform/sform if these are not full rank
% then something bad has probably happened--this will most
% likely later cause mrAlign to choke--let the user know
% here that something has gone wrong. For now, allow the user
% to continue (since they may want to manually set alignment
if rank(hdr.qform44) ~= 4
  mrWarnDlg('(mrAlignGUI) Volume qform is not full rank (This is probably a corrupted file');
end
if rank(hdr.sform44) ~= 4
  mrWarnDlg('(mrAlignGUI) Volume sform is not full rank (This is probably a corrupted file');
end

% Handle 4D file
volumeDimension = length(size(vData));
if (volumeDimension == 4)
    paramsInfo = {{'frameNum',0,'incdec=[-1 1]',sprintf('minmax=[0 %i]',hdr.dim(5)),'This volume is a 4D file, to display it as an anatomy you need to choose a particular time point or take the mean over all time points. Setting this value to 0 will compute the mean, otherwise you can select a particular timepoint to display'}};
    params = mrParamsDialog(paramsInfo,'Choose which frame of 4D file. 0 for mean');
    if isempty(params)
      return
    end
    if params.frameNum == 0
      vData = nanmean(vData,4);
    else
      vData = vData(:,:,:,params.frameNum);
    end
end

% Warning if no (qform) alignment information in the header.
% qform is initialized to identity by default in mlrImageReadNiftiHeader.
if ~(hdr.qform_code)
    mrWarnDlg('(mrAlignGUI) No alignment information in the volume header.');
end

% Warning if no (sform) base coordinate frame in the header.
% sform is initialized to identity by default in mlrImageReadNiftiHeader.
if ~(hdr.sform_code)
    mrWarnDlg('(mrAlignGUI) No base coordinate frame (i.e. sform_code = 0) in the volume header. Usually this is because the volume is straight off the magnet and has never had its base coordinate frame set. If this is a canonical base volume which you wish to align inplanes to, you should Set Base Coordinate Frame from the File menu.');
end

% Extract permutation matrix to keep track of slice orientation
permutationMatrix = getPermutationMatrix(hdr);

% Update ALIGN structure and GUI
ALIGN.volume = vData;
ALIGN.volumeHdr = hdr;
ALIGN.volumePermutation = permutationMatrix;
ALIGN.volSize = size(ALIGN.volume);
ALIGN.volumeVoxelSize = hdr.pixdim([2,3,4]);
ALIGN.coords = min(ALIGN.coords,ALIGN.volSize);
ALIGN.volumeClip = clipRange(ALIGN.volume);

% If both inplane and volume are loaded, then use the sforms from each for
% the alignment. Otherwise, use identity.
ALIGN.xform = eye(4);
if ~isempty(ALIGN.volumeHdr) & ~isempty(ALIGN.inplaneHdr)
   if ALIGN.inplaneHdr.sform_code && hdr.sform_code
      ALIGN.xform = ALIGN.volumeHdr.sform44 \ ALIGN.inplaneHdr.sform44;
   elseif ALIGN.inplaneHdr.qform_code && hdr.qform_code
      ALIGN.xform = ALIGN.volumeHdr.qform44 \ ALIGN.inplaneHdr.qform44;
   end
end

% Refresh GUI
setAlignGUI(handles,'rot',[0 0 0]);
setAlignGUI(handles,'trans',[0 0 0]);
ALIGN.guiXform = getGuiXform(handles);
setAlignGUI(handles,'nSlices',ALIGN.volSize(ALIGN.sliceOrientation));
set(handles.sliceSlider,'value',ALIGN.coords(ALIGN.sliceOrientation));
sagittalRadioButton_Callback(hObject, eventdata, handles);
refreshAlignDisplay(handles);
setAlignTitle(handles);

% Turn on Talairach option
set(handles.setTalXform,'Enable','on');
% and ability to set base coordinate frame
set(handles.setBaseCoordinateFrameMenuItem,'Enable','on');

function cRange = clipRange(image)
% Choose clipping based on histogram
histThresh = length(image(:))/1000;
[cnt, val] = hist(image(:),100);
goodVals = find(cnt>histThresh);
clipMin = val(min(goodVals));
clipMax = val(max(goodVals));
% Handle degenerate case such as image filled with NaNs
if isempty(clipMin)
    clipMin = 0;
    clipMax = 1;
end
% Handle degenerate case in which clipMax = clipMin
if (clipMax == clipMin)
    clipMax = clipMax+1;
end
cRange = [clipMin,clipMax];


% --------------------------------------------------------------------
function loadInplaneMenuItem_Callback(hObject, eventdata, handles)
global ALIGN

% Prompt user to choose inplanes. 
initPath = pwd;
pathStr = mlrGetPathStrDialog(initPath,'Choose inplane anatomy file',{'*.img;*.nii','NIFTI Files'});
if isempty(pathStr),return,end

% Load inplane file and header
h = mrMsgBox('Loading inplanes. Please wait');
[vData,hdr] = mlrImageReadNifti(pathStr);
mrCloseDlg(h);
if isempty(vData),return,end

if size(vData,3)==1
  mrWarnDlg('(mrAlignGUI) Single-slice images are not supported');
  return
end

ALIGN.inplanePath = pathStr;
% Load associated .mat file for the inplanes (used in mrLR for bases,
% and needed here for keeping track of Talairach information). If
% no such .mat file exists, create it. (If user doesn't use mrLR,
% they can just ignore this extra file
matFilename = sprintf('%s.mat',stripext(pathStr));
if ~(exist(matFilename,'file')==2) % if base doesn't already have an associated mat file
  ALIGN.inplaneBase.name = pathStr;
  ALIGN.inplaneBase.data = vData;
  ALIGN.inplaneBase.hdr = hdr;
  ALIGN.inplaneBase.permutationMatrix = getPermutationMatrix(hdr);
  ALIGN.inplaneBase.vol2mag = []; % will inherit these from the destination volume
  ALIGN.inplaneBase.vol2tal = [];
else % if there already is a base file
  load(matFilename); % then load it;
  if (exist('base')==1)  % check that this mat file actually has a base saved to it;
    ALIGN.inplaneBase = base; clear base;
  else % if some other mat file, make the base
    ALIGN.inplaneBase.name = pathStr;
    ALIGN.inplaneBase.data = vData;
    ALIGN.inplaneBase.hdr = hdr;
    ALIGN.inplaneBase.permutationMatrix = getPermutationMatrix(hdr);
    ALIGN.inplaneBase.vol2mag = []; % will inherit these from the destination volume
    ALIGN.inplaneBase.vol2tal = [];
  end
end
[tf ALIGN.inplaneBase] = isbase(ALIGN.inplaneBase); % make sure it has all the right fields;


% check rank of qform/sform if these are not full rank
% then something bad has probably happened--this will most
% likely later cause mrAlign to choke--let the user know
% here that something has gone wrong. For now, allow the user
% to continue (since they may want to manually set alignment
if rank(hdr.qform44) ~= 4
  mrWarnDlg('(mrAlignGUI) Volume qform is not full rank (This is probably a corrupted file');
end

if hdr.sform_code && (rank(hdr.sform44) ~= 4)
  mrWarnDlg('(mrAlignGUI) Volume sform is not full rank (This is probably a corrupted file');
end

% Handle 4D file
inplaneDimension = length(size(vData));
if (inplaneDimension == 4)
    paramsInfo = {{'frameNum',0,'incdec=[-1 1]',sprintf('minmax=[0 %i]',hdr.dim(5)),'This volume is a 4D file, to display it as an anatomy you need to choose a particular time point or take the mean over all time points. Setting this value to 0 will compute the mean, otherwise you can select a particular timepoint to display'}};
    params = mrParamsDialog(paramsInfo,'Choose which frame of 4D file. 0 for mean');
    if isempty(params)
      return
    end
    if params.frameNum == 0
      vData = nanmean(vData,4);
    else
      vData = vData(:,:,:,params.frameNum);
    end
end

% Warning if no (qform) alignment information in the header.
% qform is initialized to identity by default in mlrImageReadNiftiHeader.
if ~(hdr.qform_code)
    mrWarnDlg('(mrAlignGUI) No scanner alignment information in the inplane header.');
end

% Warning if no (sform) base coordinate frame in the header.
% sform is initialized to identity by default in mlrImageReadNiftiHeader.
if ~(hdr.sform_code)
    mrWarnDlg('(mrAlignGUI) No base coordinate frame (i.e. sform_code = 0) in the inplane header. Usually this is because this inplane has not yet been aligned.');
end

% Update ALIGN structure
ALIGN.inplanes = vData;
ALIGN.inplanesClip = clipRange(ALIGN.inplanes);
ALIGN.inplaneHdr = hdr;
ALIGN.inplaneSize = size(ALIGN.inplanes);
ALIGN.inplaneVoxelSize = hdr.pixdim([2,3,4]);

% If both inplane and volume are loaded, then use the sforms from each for
% the alignment. Otherwise, use identity.
ALIGN.xform = eye(4);
if ~isempty(ALIGN.volumeHdr) && ~isempty(ALIGN.inplaneHdr)
   if hdr.sform_code && ALIGN.volumeHdr.sform_code
      ALIGN.xform = ALIGN.volumeHdr.sform44 \ ALIGN.inplaneHdr.sform44;
   elseif hdr.qform_code && ALIGN.volumeHdr.qform_code
      ALIGN.xform = ALIGN.volumeHdr.qform44 \ ALIGN.inplaneHdr.qform44;
      mrWarnDlg('(mrAlignGUI) Initializing transformation from qforms.');
  else
    mrWarnDlg('(mrAlignGUI) Initializing transformation to identity.');
   end
else
  mrWarnDlg('(mrAlignGUI) Initializing transformation to identity.');
end
ALIGN.xformICCorrection = ALIGN.xform;
ALIGN.correctedInplanes = [];    
ALIGN.correctedVolume = [];      
ALIGN.xformIsRigidBody = 1; %or is it ? would be better to find out with some clever calculation...

%reset crop region
ALIGN.crop = [];

% allow source to be loaded as destination
set(handles.loadSourceAsDestination,'Enable','on');

% Refresh GUI
setAlignGUI(handles,'rot',[0 0 0]);
setAlignGUI(handles,'trans',[0 0 0]);
ALIGN.guiXform = getGuiXform(handles);
refreshAlignDisplay(handles);
setAlignTitle(handles);

% --------------------------------------------------------------------
function loadSourceAsDestination_Callback(hObject, eventdata, handles)
% hObject    handle to loadSourceAsDestination (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global ALIGN
mrAlignLoadVol(ALIGN.inplanePath,hObject,eventdata,handles);

% --------------------------------------------------------------------
function saveAlignMenuItem_Callback(hObject, eventdata, handles)
global ALIGN

if ~isfield(ALIGN,'inplaneHdr') || isempty(ALIGN.inplaneHdr)
  mrWarnDlg('(mrAlign) You have not yet loaded an inplane - no alignement to save');
  return
else
  disp('(mrAlign) Saving alignment to file...');
  sform = ALIGN.volumeHdr.sform44 * ALIGN.guiXform * ALIGN.xform;
  ALIGN.inplaneHdr = cbiSetNiftiSform(ALIGN.inplaneHdr,sform);
  ALIGN.inplaneHdr.sform_code = ALIGN.volumeHdr.sform_code;
  hdr = mlrImageWriteNiftiHeader(ALIGN.inplaneHdr,ALIGN.inplanePath);

  %also save the base structure with the appropriate vol2mag and vol2tal
  ALIGN.inplaneBase.vol2mag = ALIGN.volBase.vol2mag; % inherit from the volume
  ALIGN.inplaneBase.vol2tal = ALIGN.volBase.vol2tal; % inherit from the volume
  base = ALIGN.inplaneBase;
  matFilename = sprintf('%s.mat',stripext(base.name));
  base.data = [];base.hdr = [];
  inplanePath=fileparts(ALIGN.inplanePath);
  eval(sprintf('save %s base',fullfile(inplanePath,matFilename)));
  clear base
end
% --------------------------------------------------------------------
function saveAlignToFileMenuItem_Callback(hObject, eventdata, handles)
global ALIGN

% Extract sform
sform = ALIGN.volumeHdr.sform44 * ALIGN.guiXform * ALIGN.xform;
sform_code = ALIGN.volumeHdr.sform_code;

% Prompt user for filename(s)
pathStr = mlrGetPathStrDialog(pwd,'Choose one or more nifti files',{'*.img;*.nii','NIFTI Files'},'on');
if ~iscell(pathStr)
	pathStr = {pathStr};
end

% Loop through files and add sform to the headers
for p = 1:length(pathStr)
	if exist(pathStr{p},'file')
		hdr = mlrImageReadNiftiHeader(pathStr{p});
		hdr = cbiSetNiftiSform(hdr,sform);
                hdr.sform_code = sform_code;
		hdr = mlrImageWriteNiftiHeader(hdr,pathStr{p});
	else
		mrWarnDlg(['File ',pathStr{p},' not found.']);
	end
end

% --------------------------------------------------------------------
function saveAlignedSourceMenuItem_Callback(hObject, eventdata, handles)
global ALIGN
if isempty(ALIGN.inplanes)
  mrErrorDlg('saveAlignedSource: Must load source (inplanes) first before it can be resampled and saved');
  return
end

% Prompt user for filename(s)
pathStr = putPathStrDialog(pwd,'Specify a nifti filename',{'*.img;*.nii','NIFTI Files'});
if isempty(pathStr)
  mrWarnDlg('(mrAlignGUI) saveAlignedSource aborted');
  return
end

% Reload
[data,hdr] = mlrImageReadNifti(ALIGN.inplanePath);

% Compose xform
xform = ALIGN.guiXform * ALIGN.xform;

% Handle 3D vs 4D file
dimension = length(size(data));
volSize = ALIGN.volSize;
switch dimension
  case 3
    % Call interpVolume, but using the inverse of 'xform' (compare
    % 'transformInplanes' with 'interpVolume.m'
    interpData = interpVolume(data,inv(xform),volSize);
  case 4
    nFrames = size(data,4);
    interpData = zeros([volSize,nFrames]);
    wbh = mrWaitBar(0,'Computing alignment...');
    for f = 1:size(data,4);
      mrWaitBar(f/nFrames,wbh);
      interpData(:,:,:,f) = interpVolume(data(:,:,:,f),inv(xform),volSize);
    end
    mrCloseDlg(wbh);
  otherwise
    mrErrorDlg('saveAlignedSource: Invalid source volume, must be either 3D or 4D');
end        

% Change header
hdr = cbiSetNiftiSform(hdr,ALIGN.volumeHdr.sform44);
hdr = cbiSetNiftiQform(hdr,ALIGN.volumeHdr.qform44);

% Save and clear
mlrImageWriteNifti(pathStr,interpData,hdr);
clear interpData data

% also save the base structure with the appropriate vol2mag and vol2tal
%also save the base structure with the appropriate vol2mag and vol2tal
ALIGN.inplaneBase.vol2mag = ALIGN.volBase.vol2mag; % inherit from the volume
ALIGN.inplaneBase.vol2tal = ALIGN.volBase.vol2tal; % inherit from the volume
base = ALIGN.inplaneBase;
matFilename = sprintf('%s.mat',stripext(base.name));
base.data = [];base.hdr = [];
eval(sprintf('save %s base',matFilename))
clear base
% --------------------------------------------------------------------
function saveResampledDestination_Callback(hObject, eventdata, handles)
% hObject    handle to saveResampledDestination (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global ALIGN
if isempty(ALIGN.inplanes)
  mrErrorDlg('saveAlignedSource: Must load source (inplanes) first before it can be resampled and saved');
  return
end

% Prompt user for filename(s)
pathStr = putPathStrDialog(pwd,'Specify a nifti filename',{'*.img;*.nii','NIFTI Files'});
if isempty(pathStr)
  mrWarnDlg('(mrAlignGUI) saveAlignedSource aborted');
  return
end

% set output voxel size to be the same as the volume anatomy from which we are sampling
outputVoxSize = [ALIGN.volumeVoxelSize(1) ALIGN.volumeVoxelSize(2) ALIGN.inplaneVoxelSize(3)];

% get transformation for inplane voxel size to volume voxel size, this creates images with
% voxel sizes the have the same inplane dimension as the volume anatomy, but with slice
% thickness the same as the inplanes.
voxSizeXform = diag([outputVoxSize(:)'./ALIGN.inplaneVoxelSize(:)' 1]);

% Compose xform
xform = ALIGN.guiXform * ALIGN.xform * voxSizeXform;

% get the size of the inplanes fov in the inplane dimensions
inplaneFOV = ALIGN.inplaneSize(1:2).*ALIGN.inplaneVoxelSize(1:2)';

% compute the volume size we will need to cover the same area
volSize = round([inplaneFOV./outputVoxSize(1:2)]);

% In the slice dimension we keep the same volume size as the inplanes
volSize(3) = ALIGN.inplaneSize(3);

% Call interpVolume
interpData = interpVolume(ALIGN.volume,xform,volSize);

% now create the nifti header
hdr = cbiCreateNiftiHeader(interpData);
hdr = cbiSetNiftiQform(hdr,ALIGN.volumeHdr.qform44*xform);
hdr = cbiSetNiftiSform(hdr,ALIGN.volumeHdr.qform44*xform);
hdr.pixdim(2:4) = outputVoxSize;

% Save and clear
[byteswritten,hdr] = mlrImageWriteNifti(pathStr,interpData,hdr);
clear interpData data

% also save the base structure with the appropriate vol2mag and vol2tal
%also save the base structure with the appropriate vol2mag and vol2tal
ALIGN.inplaneBase.vol2mag = ALIGN.volBase.vol2mag; % inherit from the volume
ALIGN.inplaneBase.vol2tal = ALIGN.volBase.vol2tal; % inherit from the volume
base = ALIGN.inplaneBase;
matFilename = sprintf('%s.mat',stripext(base.name));
base.data = [];base.hdr = [];
eval(sprintf('save %s base',matFilename))
clear base

% --------------------------------------------------------------------
function importMenu_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
% Imports xform from separate file
function importAlignMenuItem_Callback(hObject, eventdata, handles)
global ALIGN

% Prompt user for file and load it
pathstr = mlrGetPathStrDialog(pwd,'Choose alignment file','*.mat');
if ~exist(pathstr,'file')
    return
end
load(pathstr);
% Warning if old version.
if ~exist('mrAlignVersion','var') | ~isnumeric(mrAlignVersion) | (mrAlignVersion < ALIGN.version)
    mrWarnDlg('(mrAlignGUI) This alignment file appears to correspond to an older version of mrAlign. You may need to use "Import" from the "File" menu.');
end
% Error if xform isn't loaded from the file.
if ~exist('xform','var')
    mrErrorDlg('Invalid alignment file.');
end

% Update ALIGN structure and GUI.
ALIGN.xform = xform;
setAlignGUI(handles,'rot',[0 0 0]);
setAlignGUI(handles,'trans',[0 0 0]);
ALIGN.guiXform = getGuiXform(handles);
refreshAlignDisplay(handles);

% --------------------------------------------------------------------
% Imports xform from separate file
function importAlignFromFlirtMenuItem_Callback(hObject, eventdata, handles)
global ALIGN

if isempty(ALIGN.volumeVoxelSize)
  mrWarnDlg('(mrAlignGUI) You must load source and destination volumes before importing an alignment');
  return
end

% Prompt user for file and load it
pathstr = mlrGetPathStrDialog(pwd,'Choose FLIRT alignment text file','*.*');
if ~exist(pathstr,'file')
    return
end

%we assume that the file is a 4*4 float matrix in an ascii file
try
  xform=dlmread(pathstr);
catch id
    mrWarnDlg('(mrAlignGUI) Invalid FLIRT alignment file.');
    disp(id.message);
    return;
end

% Check that file contains a 4*4 matrix with a 1 in the last cell
if ~all(size(xform)==4) || abs(xform(4,4)-1)>1e-7
    mrWarnDlg('(mrAlignGUI) Invalid transformation matrix in FLIRT alignment file.');
    disp(xform);
    return;
end

%translation values are in mm and must be converted to voxels of the destination image !!!
xform(1:3,4) = xform(1:3,4)./ALIGN.volumeVoxelSize;

% Update ALIGN structure and GUI.
ALIGN.xform = xform;
setAlignGUI(handles,'rot',[0 0 0]);
setAlignGUI(handles,'trans',[0 0 0]);
ALIGN.guiXform = getGuiXform(handles);
refreshAlignDisplay(handles);

% --------------------------------------------------------------------
function composeAlignmentMenuItem_Callback(hObject, eventdata, handles)
global ALIGN

% Prompt user for first alignment file and load it.
pathstr = mlrGetPathStrDialog(pwd,'Choose first alignment file','*.mat');
if ~exist(pathstr,'file')
    return
end
load(pathstr);
% Warning if old version.
if ~exist('mrAlignVersion','var') | ~isnumeric(mrAlignVersion) | (mrAlignVersion < ALIGN.version)
    mrWarnDlg('(mrAlignGUI) This alignment file appears to correspond to an older version of mrAlign.');
end
% Error if  xform isn't loaded from the file.
if ~exist('xform','var')
    mrErrorDlg('Invalid alignment file.');
else
    xform1 = xform;
end

% Prompt user for second alignment file and load it.
pathstr = mlrGetPathStrDialog(pwd,'Choose second alignment file','*.mat');
load(pathstr);
if ~exist(pathstr,'file')
    return
end
% Warning if old version.
if ~exist('mrAlignVersion','var') | ~isnumeric(mrAlignVersion) | (mrAlignVersion < ALIGN.version)
    mrWarnDlg('(mrAlignGUI) This alignment file appears to correspond to an older version of mrAlign.');
end  
% Error if xform isn't loaded from the file.
if ~exist('xform','var')
    mrErrorDlg('Invalid alignment file.');
else
    xform2 = xform;
end

% Compose the alignments. Update ALIGN structure and GUI.
ALIGN.xform = xform2 * xform1;
setAlignGUI(handles,'rot',[0 0 0]);
setAlignGUI(handles,'trans',[0 0 0]);
ALIGN.guiXform = getGuiXform(handles);
refreshAlignDisplay(handles);

% --------------------------------------------------------------------
% Imports alignment from mrAlign-4.2 or earlier
function importOldAlignMenuItem_Callback(hObject, eventdata, handles)
global ALIGN

% Prompt user and load it
pathstr = mlrGetPathStrDialog(pwd,'Choose alignment file','*.mat');
if ~exist(pathstr,'file')
    return
end
load(pathstr);
% Error if  rot, trans, and scaleFac aren't loaded from the file.
if ~exist('rot','var') | ~exist('trans','var') | ~exist('scaleFac','var')
    myErrorDlg('Invalid alignment file.');
end

% Use loaded rot, trans, and scaleFac to compute xform.
S1 = [diag(scaleFac(1,:)) zeros(3,1); 0 0 0 1];
S2 = [diag(scaleFac(2,:)) zeros(3,1); 0 0 0 1];
Mi = [rot trans(:); 0 0 0 1];
xform = S2*Mi*inv(S1);

% Convert from mrAlign-4.2 convention to current convention (permute the
% dimensions and flip two of them).
% *** This hasn't been exhaustively test. Possible bug is that volSize(2)
% and volSize(3) may need to be swapped below.
xform = [xform(3,:); xform(1,:); xform(2,:); xform(4,:)];
xform(2,:) = -xform(2,:);
xform(2,4) = xform(2,4) + ALIGN.volSize(2);
xform(3,:) = -xform(3,:);
xform(3,4) = xform(3,4) + ALIGN.volSize(3);

% Shift xform: matlab indexes from 1 but nifti uses 0,0,0 as the origin. 
shiftXform = shiftOriginXform;
xform = inv(shiftXform)* xform;

% Update ALIGN structure and GUI
ALIGN.xform = xform;
setAlignGUI(handles,'rot',[0 0 0]);
setAlignGUI(handles,'trans',[0 0 0]);
ALIGN.guiXform = getGuiXform(handles);
refreshAlignDisplay(handles);

% --------------------------------------------------------------------
function exportMenu_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
% Exports xform to separate file
function exportAlignMenuItem_Callback(hObject, eventdata, handles)
global ALIGN

xform = ALIGN.guiXform * ALIGN.xform;
mrAlignVersion = ALIGN.version;
pathstr = putPathStrDialog(pwd,'Specify alignment file','*.mat');
% pathstr = [] if aborted
if ~isempty(pathstr)
    save(pathstr,'xform','mrAlignVersion');
end

% --------------------------------------------------------------------
% Exports alignment for mrLoadRet-3.1 or earlier 
function exportOldAlignMenuItem_Callback(hObject, eventdata, handles)
global ALIGN

% Extract xform
xform = ALIGN.guiXform * ALIGN.xform;

% Shift xform: matlab indexes from 1 but nifti uses 0,0,0 as the origin. 
shiftXform = shiftOriginXform;
xform = xform * shiftXform;

% Reverse of the logic in importAlignMenuItem_Callback (mrAlignGUI) to
% convert from the current convention back to mrAlign-4.2 convention.
% *** This hasn't been exhaustively test. Possible bug is that volSize(2)
% and volSize(3) may need to be swapped below.
xform(3,4) = xform(3,4) - ALIGN.volSize(3);
xform(3,:) = -xform(3,:);
xform(2,4) = xform(2,4) - ALIGN.volSize(2);
xform(2,:) = -xform(2,:);
xform = [xform(2,:); xform(3,:); xform(1,:); xform(4,:)];

% Get scaleFac from voxel sizes
% *** This hasn't been exhaustively test. Possible bug is that x and y
% voxel sizes need to be swapped.
scaleFac = [1./ALIGN.inplaneVoxelSize'; 1./ALIGN.volumeVoxelSize'];

% compute rot and trans from 4x4
b = (xform(1:3,4))';
A = xform(1:3,1:3);
trans = b ./ scaleFac(2,:);
rot = diag(scaleFac(2,:)) \ A / diag( 1./scaleFac(1,:));

% Label it with version number. 
mrAlignVersion = ['exported from mrAlign ',num2str(ALIGN.version)];

pathstr = putPathStrDialog(pwd,'Specify alignment file','*.mat');
% pathstr = [] if aborted
if ~isempty(pathstr)
    save(pathstr,'rot','trans','scaleFac','mrAlignVersion');
end

% --------------------------------------------------------------------
function setBaseCoordinateFrameMenuItem_Callback(hObject, eventdata, handles)
global ALIGN

% check sform_code
sform_code = ALIGN.volumeHdr.sform_code;

if sform_code == 3 % just set the vol2tal and vol2mag fields
  disp(sprintf('(mrAlignGUI) Talairach volume (sform_code=3). Setting vol2tal and vol2mag fields'));
  ALIGN.volBase.vol2tal = ALIGN.volumeHdr.sform44;
  ALIGN.volBase.vol2mag = ALIGN.volumeHdr.qform44;
else % if sform_code is 0 or 1, set sform = qform and get vol2mag.
  sform = ALIGN.volumeHdr.qform44;
  ALIGN.volumeHdr = cbiSetNiftiSform(ALIGN.volumeHdr,sform);
  hdr = mlrImageWriteNiftiHeader(ALIGN.volumeHdr,ALIGN.volumePath);
  ALIGN.volBase.vol2mag = sform; % set the vol2mag field
  disp(sprintf('(mrAlignGUI) Setting sform_code=1 and sform equal to qform and saving in nifti header'));
end
  
%also save the matFile for the base, to save vol2mag and vol2tal
base = ALIGN.volBase;
matFilename = sprintf('%s.mat',stripext(base.name));
base.data = [];base.hdr = [];
disp(sprintf('(mrAlignGUI) Saving %s file',matFilename));
volumePath=fileparts(ALIGN.volumePath);
eval(sprintf('save %s base',fullfile(volumePath,matFilename)));
clear base % so don't get confused with the inplane base


% --------------------------------------------------------------------
function writeTifMenuItem_Callback(hObject, eventdata, handles)
pathstr = putPathStrDialog(pwd,'Specify alignment file','*.tif');
% pathstr = [] if aborted
if ~isempty(pathstr)
	img = refreshAlignDisplay(handles);
	imwrite(img,pathstr,'tif');
end

% --------------------------------------------------------------------
function quitMenuItem_Callback(hObject, eventdata, handles)
clear global ALIGN

% remember figure location
mrSetFigLoc('mrAlignGUI',get(handles.figure1,'Position'));
% close figure
delete(handles.figure1);
% save .mrDefaults in the home directory
saveMrDefaults;



% --------------------------------------------------------------------
function editMenu_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function prefMenu_Callback(hObject, eventdata, handles)

mrEditPrefs;

% --------------------------------------------------------------------
function manualAlignmentMenu_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function dummy = initializeIdentityMenuItem_Callback(hObject, eventdata, handles)
dummy=[]; %this is so this callback can be called through mrParamsDialog
global ALIGN

% Set xform to identity, but scaled by voxel sizes
% *** Not yet test/debugged ***

%first check the qform matrices of the inplane and volumes (if loaded)
if isempty(ALIGN.volumeHdr) || isempty(ALIGN.inplaneHdr)
	mrWarnDlg('(mrAlignGUI) Load source and destination');
	return
else
   answer = '';
   if ALIGN.inplaneHdr.sform_code && ALIGN.volumeHdr.sform_code && det(ALIGN.inplaneHdr.sform44)*det(ALIGN.volumeHdr.sform44)<0
      question = 'It seems that one of the axes of this source has been previously flipped in order for the data to be in right-handed space.';
      question = [question ' Setting the transformation matrix to identity will undo this.']; 
      question = [question ' You should flip one axis before continuing.']; 
      question = [question ' What do you want to do ?']; 
      answer = questdlg(question,'Space Orientation Mismatch ?','Flip X', 'Do nothing', 'Do nothing');
   end
end

% jg: These were set to 1/voxelSize which I believe
% was wrong. Flipped these and now this seems to work
% correctly
inplaneXform = eye(4);
inplaneXform(1,1) = ALIGN.inplaneVoxelSize(1);
inplaneXform(2,2) = ALIGN.inplaneVoxelSize(2);
inplaneXform(3,3) = ALIGN.inplaneVoxelSize(3);
volumeXform = eye(4);
volumeXform(1,1) = ALIGN.volumeVoxelSize(1);
volumeXform(2,2) = ALIGN.volumeVoxelSize(2);
volumeXform(3,3) = ALIGN.volumeVoxelSize(3);
ALIGN.xform = volumeXform\inplaneXform;
ALIGN.xformIsRigidBody = 1;

switch answer
   case 'Flip X'
      flipXMenuItem_Callback([], [], handles);
   case 'Flip Y'
      flipYMenuItem_Callback([], [], handles);
   case 'Flip Z'
      flipZMenuItem_Callback([], [], handles);
end    


% Reset GUI
setAlignGUI(handles,'rot',[0 0 0]);
setAlignGUI(handles,'trans',[0 0 0]);
ALIGN.guiXform = getGuiXform(handles);
refreshAlignDisplay(handles);

% --------------------------------------------------------------------
function flipXMenuItem_Callback(hObject, eventdata, handles)
global ALIGN
xform = ALIGN.guiXform * ALIGN.xform;
ALIGN.xform = xform * [-1 0 0 ALIGN.inplaneSize(1); 0 1 0 0; 0 0 1 0; 0 0 0 1];
setAlignGUI(handles,'rot',[0 0 0]);
setAlignGUI(handles,'trans',[0 0 0]);
ALIGN.guiXform = getGuiXform(handles);
refreshAlignDisplay(handles);

% --------------------------------------------------------------------
function flipYMenuItem_Callback(hObject, eventdata, handles)
global ALIGN
xform = ALIGN.guiXform * ALIGN.xform;
ALIGN.xform = xform * [1 0 0 0; 0 -1 0 ALIGN.inplaneSize(2); 0 0 1 0; 0 0 0 1];
setAlignGUI(handles,'rot',[0 0 0]);
setAlignGUI(handles,'trans',[0 0 0]);
ALIGN.guiXform = getGuiXform(handles);
refreshAlignDisplay(handles);

% --------------------------------------------------------------------
function flipZMenuItem_Callback(hObject, eventdata, handles)
global ALIGN
xform = ALIGN.guiXform * ALIGN.xform;
ALIGN.xform = xform * [1 0 0 0; 0 1 0 0; 0 0 -1 ALIGN.inplaneSize(3); 0 0 0 1];
setAlignGUI(handles,'rot',[0 0 0]);
setAlignGUI(handles,'trans',[0 0 0]);
ALIGN.guiXform = getGuiXform(handles);
refreshAlignDisplay(handles);


% --------------------------------------------------------------------
function swapXYmenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to swapXYmenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global ALIGN
xform = ALIGN.guiXform * ALIGN.xform;
ALIGN.xform = xform * [0 1 0 0; 1 0 0 0; 0 0 1 0; 0 0 0 1];
setAlignGUI(handles,'rot',[0 0 0]);
setAlignGUI(handles,'trans',[0 0 0]);
ALIGN.guiXform = getGuiXform(handles);
refreshAlignDisplay(handles);


% --------------------------------------------------------------------
function swapYZmenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to swapYZmenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global ALIGN
xform = ALIGN.guiXform * ALIGN.xform;
ALIGN.xform = xform * [1 0 0 0; 0 0 1 0; 0 1 0 0; 0 0 0 1];
setAlignGUI(handles,'rot',[0 0 0]);
setAlignGUI(handles,'trans',[0 0 0]);
ALIGN.guiXform = getGuiXform(handles);
refreshAlignDisplay(handles);


% --------------------------------------------------------------------
function swapXZmenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to swapXZmenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global ALIGN
xform = ALIGN.guiXform * ALIGN.xform;
ALIGN.xform = xform * [0 0 1 0; 0 1 0 0; 1 0 0 0; 0 0 0 1];
setAlignGUI(handles,'rot',[0 0 0]);
setAlignGUI(handles,'trans',[0 0 0]);
ALIGN.guiXform = getGuiXform(handles);
refreshAlignDisplay(handles);
% --------------------------------------------------------------------
function computeAlignmentMenu_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function initializeFromQformMenuItem_Callback(hObject, eventdata, handles)
global ALIGN

if ~isfield(ALIGN.volumeHdr,'qform44') || ~isfield(ALIGN.inplaneHdr,'qform44')
  mrWarnDlg('(mrAlignGUI) Need to load both src and dest images');
  return
end

% Error if there's no alignment information in the header.
% This would happen if these were analyze, not nifti, files.
if isempty(ALIGN.volumeHdr.qform44)
    mrErrorDlg('No alignment information in the volume header.');
end
if ~isempty(ALIGN.volumeHdr.qform_code) && ~ALIGN.volumeHdr.qform_code
    mrWarnDlg('(mrAlignGUI) Volume qform_code is not set.');
end
if isempty(ALIGN.inplaneHdr.qform44)
    mrErrorDlg('No alignment information in the inplane header.');
end
if ~isempty(ALIGN.inplaneHdr.qform_code) && ~ALIGN.inplaneHdr.qform_code
    mrWarnDlg('(mrAlignGUI) Inplanes qform_code is not set.');
end

% Compute alignment by composing the qforms from the two nifti headers.
ALIGN.xform = ALIGN.volumeHdr.qform44 \ ALIGN.inplaneHdr.qform44;

% Reset GUI
setAlignGUI(handles,'rot',[0 0 0]);
setAlignGUI(handles,'trans',[0 0 0]);
ALIGN.guiXform = getGuiXform(handles);
refreshAlignDisplay(handles);

% --------------------------------------------------------------------
function initializeFromSformMenuItem_Callback(hObject, eventdata, handles)
global ALIGN

if ~isfield(ALIGN.volumeHdr,'sform44') || ~isfield(ALIGN.inplaneHdr,'sform44')
  mrWarnDlg('(mrAlignGUI) Need to load both src and dest images');
  return
end

% Error if there's no alignment information in the header.
% This would happen if these were analyze, not nifti, files.
if isempty(ALIGN.volumeHdr.sform44)
    mrErrorDlg('No alignment information in the volume header.');
end
if ~isempty(ALIGN.volumeHdr.sform_code) && ~ALIGN.volumeHdr.sform_code
    mrWarnDlg('(mrAlignGUI) Volume sform_code is not set.');
end
if isempty(ALIGN.inplaneHdr.sform44)
    mrErrorDlg('No alignment information in the inplane header.');
end
if ~isempty(ALIGN.inplaneHdr.sform_code) && ~ALIGN.inplaneHdr.sform_code
    mrWarnDlg('(mrAlignGUI) Inplanes sform_code is not set.');
end

% Compute alignment by composing the sforms from the two nifti headers.
ALIGN.xform = ALIGN.volumeHdr.sform44 \ ALIGN.inplaneHdr.sform44;

% Reset GUI
setAlignGUI(handles,'rot',[0 0 0]);
setAlignGUI(handles,'trans',[0 0 0]);
ALIGN.guiXform = getGuiXform(handles);
refreshAlignDisplay(handles);

% --------------------------------------------------------------------
function cropInplanesMenuItem_Callback(hObject, eventdata, handles)
global ALIGN
ALIGN.crop = selectCropRegion(ALIGN.inplanes);

% --------------------------------------------------------------------
function resetCropInplanesMenuItem_Callback(hObject, eventdata, handles)
global ALIGN
ALIGN.crop = [];

% --------------------------------------------------------------------
function advancedAlignmentMenuItem_Callback(hObject, eventdata, handles)
global ALIGN

paramsInfo = {...
  {'identity', 0,'type=pushbutton','buttonString=Set Alignment to Identity','callback',{@initializeIdentityMenuItem_Callback,[],[],handles},'passCallbackOutput=0',...
      'Sets the transformation to dentiy, taking into account differences in voxel sizes between the source and destination volumes.'},...
  {'qform', 0,'type=pushbutton','buttonString=Initialize Alignment from Header Qform','callback',{@initializeFromQformMenuItem_Callback,[],[],handles},'passCallbackOutput=0',...
      'Sets the transformation to the qform matrix of the NIFTI header. This will revert back to the original alignment of the volumes right off the scanner.'},...
  {'sform', 0,'type=pushbutton','buttonString=Initialize Alignment from Header Sform','callback',{@initializeFromSformMenuItem_Callback,[],[],handles},'passCallbackOutput=0',...
      'Sets the transformation to the sform matrix of the NIFTI header. This will revert back to the last saved alignment.'},...
  {'set', 0,'type=pushbutton','buttonString=Set Crop Region','callback',@setCropRegion,'passCallbackOutput=0',...
      'Sets a 3D box in the source volume that constrains the alignment estimation. Voxels outside this box are ignored in the motion estimation'},...
  {'reset', 0,'type=pushbutton','buttonString=Reset Crop Region','callback',{@resetCropInplanesMenuItem_Callback},'passCallbackOutput=0',...
      'Resets the 3D box to the whold source volume.'},...
  {'robust',ALIGN.robustAlignment,'type=checkbox','buttonString=Robust Alignment','callback',@setAlignmentOption,'callbackArg','robust',...
      'Uses robust non-linear iterative motion estimation algorithm.'},...
  {'rigid',ALIGN.rigidBodyAlignment,'type=checkbox','buttonString=Rigid-body Alignment','','callback',@setAlignmentOption,'callbackArg','rigid',...
      'Restricts the search space to rigid transformations (=6 degrees of freedom). If unchecked, an additional parameter of scaling is estiamting (7DOF)'},...
  {'reverse',ALIGN.reverseContrastAlignment,'type=checkbox','buttonString=Reverse Contrast (T2*)','callback',@setAlignmentOption,'callbackArg','reverse',...
      'Reverses the source image contrast and treats the voxel values as periodic (phases of complex number on the unit circle). This is used mainly to align high-resolution T2(*)-weighted images to T1-weighted images'},...
  {'ignore',ALIGN.ignoreZeroVoxels,'type=checkbox','buttonString=Ignore Zero Voxels','callback',@setAlignmentOption,'callbackArg','ignore',...
      'Ignores voxels that are exactly 0 in the source and destination volumes'},...
  {'display',ALIGN.displayAlignmentSteps,'type=checkbox','buttonString=Display Alignment Steps','callback',@setAlignmentOption,'callbackArg','display',...
      'Displays the new alignment in the main window at each iteration of the computation'},...
  {'coarse', 0,'type=pushbutton','buttonString=Compute Coarse Alignment','callback',{@computeAlignment,'coarse',handles},'passCallbackOutput=0',...
      'Downsamples the images to a minimum of 3 mm resolution in all dimensions (2mm for reverse contrast) and computes the alignment'},...
  {'fine', 0,'type=pushbutton','buttonString=Compute Fine Alignment','callback',{@computeAlignment,'fine',handles},'passCallbackOutput=0',...
      'Computes the alignment at the native volume resolution (to a 1mm minimum voxels; .5 for reverse contrast)'},...
  {'stop', 0,'type=pushbutton','buttonString=Stop Computing','callback',{@stopComputingMenuItem_Callback},'passCallbackOutput=0',...
      'Stops the current computation'},...
  {'undo', 0,'type=pushbutton','buttonString=Undo Last Alignment','callback',{@undoLastAlignmentMenuItem_Callback,[],[],handles},'passCallbackOutput=0',...
      'Undoes the last computed alignment and restore any manual alignment pre-existing to this alignment computation.'},...
   };

  % display dialog
  mrParamsDialog(paramsInfo,'Advanced Alignment Menu','modal=0');
  
function setCropRegion
  %this function is just so that we get the new inplane volume if it has changed
  global ALIGN
  ALIGN.crop=selectCropRegion(ALIGN.inplanes);
  
function setAlignmentOption(option,params)
global ALIGN

switch option
  case 'robust'
    ALIGN.robustAlignment = params.robust;
  case 'rigid'
    ALIGN.rigidBodyAlignment = params.rigid;
  case 'reverse'
    ALIGN.reverseContrastAlignment = params.reverse;
  case 'ignore'
    ALIGN.ignoreZeroVoxels = params.ignore;
  case 'display'
    ALIGN.displayAlignmentSteps = params.display;
end    

% --------------------------------------------------------------------
function coarseAlignmentMenuItem_Callback(hObject, eventdata, handles)

computeAlignment('coarse',handles);

% --------------------------------------------------------------------
function fineAlignmentMenuItem_Callback(hObject, eventdata, handles)

computeAlignment('fine',handles);

% --------------------------------------------------------------------
function stopComputingMenuItem_Callback(hObject, eventdata, handles)

global ALIGN
if ALIGN.currentlyComputingAlignment
   ALIGN.stopComputingAlignment = 1;
   mrWarnDlg('(mrAlignGUI) Cancelling alignment computation. This might take some time. Please wait...');
else
   mrWarnDlg('(mrAlignGUI) No Alignment Computation is running');
end

% --------------------------------------------------------------------
function undoLastAlignmentMenuItem_Callback(hObject, eventdata, handles)

global ALIGN
if isempty(ALIGN.oldXform)
  mrWarnDlg('(mrAlignGUI) There is no stored alignment')
else
  ALIGN.xform = ALIGN.oldXform;
  ALIGN.oldXform = [];
  setAlignGUI(handles,'rot',ALIGN.oldGuiRotation);
  setAlignGUI(handles,'trans',ALIGN.oldGuiTranslation);
  ALIGN.guiXform = getGuiXform(handles);
  ALIGN.oldGuiRotation = [];       
  ALIGN.oldGuiTranlastion = [];    
  refreshAlignDisplay(handles);
end

% --------------------------------------------------------------------
function mutualInformationMenuItem_Callback(hObject, eventdata, handles)

mrWarnDlg('(mrAlignGUI) Mutual information registration not yet implemented.');
return;

global ALIGN
% Options
% opts = optimset('fminunc');
opts = optimset('MaxFunEvals',100,'TolFun',1e-2,'DiffMaxChange',1,'DiffMinChange',0.1);
% Initial values (rot in deg then trans in voxels)
% *** Should be initialized to current xform, extracted from:
%     ALIGN.xform * ALIGN.guiXform;
x0 = zeros(1,6);
% Search
x = fminunc(@mutualInformationFun,x0,opts);
ALIGN.xform = extractXform(x);
ALIGN.xformICCorrection = ALIGN.xform;

% Reset GUI and refresh display
setAlignGUI(handles,'rot',[0 0 0]);
setAlignGUI(handles,'trans',[0 0 0]);
ALIGN.guiXform = getGuiXform(handles);
refreshAlignDisplay(handles);

function val = mutualInformationFun(x)
global ALIGN
display(x)
% Interpolate the volume
xform = extractXform(x);
[NyI NxI NzI] = size(ALIGN.inplanes);
interpolatedVolume = regInplanes(ALIGN.volume, NxI, NyI, NzI, xform);
% Crop
if ~isempty(ALIGN.crop)
    crop = ALIGN.crop;
    inpCrop = ALIGN.inplanes([crop(1,2):crop(2,2)],[crop(1,1):crop(2,1)],:);
    volCrop = interpolatedVolume([crop(1,2):crop(2,2)],[crop(1,1):crop(2,1)],:);
else
    inpCrop = ALIGN.crop;
    volCrop = interpolatedVolume;
end
% Compute mutual information
val = - mutualInformation(inpCrop,volCrop,1);
display(val)

function xform = extractXform(x)
global ALIGN
addXform = eye(4);
x(1:3) = x(1:3)*pi/180;
cosx = cos(x(1));		sinx = sin(x(1));
cosy = cos(x(2));		siny = sin(x(2));
cosz = cos(x(3));		sinz = sin(x(3));
addXform(1,1) = cosz*cosy+sinz*sinx*siny;
addXform(1,2) = sinz*cosy-cosz*sinx*siny;
addXform(1,3) = cosx*siny;
addXform(2,1) = -sinz*cosx;
addXform(2,2) = cosz*cosx;
addXform(2,3) = sinx;
addXform(3,1) = sinz*sinx*cosy-cosz*siny;
addXform(3,2) = -cosz*sinx*cosy-sinz*siny;
addXform(3,3) = cosx*cosy;
addXform(1,4) = x(4);
addXform(1,5) = x(5);
addXform(1,6) = x(6);
xform = addXform * ALIGN.xform;


% --------------------------------------------------------------------
function checkAlignmentMenu_Callback(hObject, eventdata, handles)


% --------------------------------------------------------------------
function showMosaicMenuItem_Callback(hObject, eventdata, handles)

global ALIGN

if isempty(ALIGN.volume) || isempty(ALIGN.inplanes)
	mrWarnDlg('(mrAlignGUI) Load Volume and Load Inplanes before checking alignment.');
	return
end
if isempty(ALIGN.xform) 
	mrWarnDlg('(mrAlignGUI) Load, Initialize, or Compute the alignment before checking it.');
	return
end

if isequal(ALIGN.xformICCorrection,ALIGN.guiXform * ALIGN.xform) &&...        % if the alignment hasn't changed (this should be the most common scenario) 
      ~isempty(ALIGN.correctedInplanes) && ~isempty(ALIGN.correctedVolume)    % AND the corrected volumes have been computed
  
   %then just make mosaic
   inplanes = ALIGN.correctedInplanes;
   volume = ALIGN.correctedVolume;
   mosaic = imageMosaic(volume,inplanes);

else  %otherwise, we have to recompute the interpolated/corrected inplanes/volume
   
   [inplanes,volume,mosaic] = checkAlignment(ALIGN.reversedContrast+1);

end
   
% open checkAlignmentGUI
checkAlignmentGUI(inplanes,mosaic,volume,ALIGN.reversedContrast+1);



% --------------------------------------------------------------------
function jointHistogramMenuItem_Callback(hObject, eventdata, handles)
global ALIGN

% Interpolate the volume
[NyI NxI NzI] = size(ALIGN.inplanes);
h = mrMsgBox('Wait while interpolating the inplanes...'); drawnow
interpolatedVolume = regInplanes(ALIGN.volume, NxI, NyI, NzI, ALIGN.xform);
mrCloseDlg(h);
% Compute histogram and display it
histogram = flipud(hist2(ALIGN.inplanes,interpolatedVolume));
FF = figure;
imagesc(log(histogram));
axis('image'); colormap('gray'); axis('off');
disp('Press a key to continue...')
pause
close(FF)


% --------------------------------------------------------------------
function exportMrLoadRet4_Callback(hObject, eventdata, handles)
% hObject    handle to exportMrLoadRet4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global ALIGN

% Extract sform
if ~isempty(ALIGN.volumeHdr)
  sform = ALIGN.volumeHdr.sform44 * ALIGN.guiXform * ALIGN.xform;
  sform_code = ALIGN.volumeHdr.sform_code;
else
  mrWarnDlg('(mrAlignGUI) Volume header does not exist yet');
  return
end
if ~isempty(ALIGN.volBase)
  vol2mag = ALIGN.volBase.vol2mag;
  vol2tal = ALIGN.volBase.vol2tal;
else
  mrWarnDlg('(mrAlignGUI) No volume loaded');
  return
end

% save sform to mrLoadRet4
saveSform(sform,sform_code,vol2tal,vol2mag);


% --------------------------------------------------------------------
function setSform_Callback(hObject, eventdata, handles)
% hObject    handle to setSform (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global ALIGN
ALIGN.guiXform = getGuiXform(handles);

paramsInfo = {};
if ~isempty(ALIGN.volumeHdr)
  paramsInfo{end+1}{1} = 'destQform';
  paramsInfo{end}{2} = ALIGN.volumeHdr.qform44;
  paramsInfo{end}{3} = 'editable=0';
  paramsInfo{end}{4} = 'Qform is the transformation to magnet coordinates';
  paramsInfo{end+1}{1} = 'destSform';
  paramsInfo{end}{2} = ALIGN.volumeHdr.sform44;
  paramsInfo{end}{3} = 'editable=0';
  paramsInfo{end}{4} = 'Sform is set by mrAlign to be to the transfomration to the coordinates of the base volume anatomy';
end
if ~isempty(ALIGN.inplaneHdr)
  paramsInfo{end+1}{1} = 'srcQform';
  paramsInfo{end}{2} = ALIGN.inplaneHdr.qform44;
  paramsInfo{end}{3} = 'editable=0';
  paramsInfo{end}{4} = 'Qform is the transformation to magnet coordinates';
  paramsInfo{end+1}{1} = 'srcSform';
  paramsInfo{end}{2} = ALIGN.volumeHdr.sform44 * ALIGN.guiXform * ALIGN.xform;

  paramsInfo{end}{3} = 'Sform is set by mrAlign to be to the transfomration to the coordinates of the base volume anatomy';
end

% put up dialog
params = mrParamsDialog(paramsInfo,'Set sform of source directly');

% user hit cancel
if isempty(params) || ~isfield(params,'srcSform')
  return
end

% set the transform 
ALIGN.inplaneHdr.sform44 = params.srcSform;
ALIGN.xform = ALIGN.volumeHdr.sform44 \ ALIGN.inplaneHdr.sform44;
setAlignGUI(handles,'rot',[0 0 0]);
setAlignGUI(handles,'trans',[0 0 0]);
ALIGN.guiXform = getGuiXform(handles);
refreshAlignDisplay(handles);


% --------------------------------------------------------------------
function setSourceToDestination_Callback(hObject, eventdata, handles)
% hObject    handle to setSourceToDestination (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


global ALIGN
ALIGN.guiXform = getGuiXform(handles);

% get srcToDestXform
paramsInfo = {};
if ~isempty(ALIGN.volumeHdr) && ~isempty(ALIGN.inplaneHdr)
  paramsInfo{end+1}{1} = 'srcToDestXform';
  paramsInfo{end}{2} = ALIGN.guiXform * ALIGN.xform;
  paramsInfo{end}{3} = 'Transformation from source to destination coordinate system';
end

% put up dialog
params = mrParamsDialog(paramsInfo,'Set xfrom of source to dest');

% user hit cancel
if isempty(params)
  return
end

% set the transform 
ALIGN.xform = params.srcToDestXform;
setAlignGUI(handles,'rot',[0 0 0]);
setAlignGUI(handles,'trans',[0 0 0]);
ALIGN.guiXform = getGuiXform(handles);
refreshAlignDisplay(handles);


% --------------------------------------------------------------------
function talMenu_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function setTalXform_Callback(hObject, eventdata, handles)
global ALIGN
% define the cannonical talairach coordinates
  tAC = [0 0 0]';
  tPC = [0 -24 0]';
  tSAC = [0 0 72]';
  tIAC = [0 0 -42]';
  tPPC = [0 -102 0]';
  tAAC = [0 68 0]';
  tLAC = [-62 0 0]';
  tRAC = [62 0 0]';
  talPoints = [tAC tPC tSAC tIAC tPPC tAAC tLAC tRAC];
  talPoints(4, :) = ones(1,size(talPoints,2));

% Initialize the 8 points used for defining the tal Transform if possible
% then run the talairach program to allow user to define the 8 points

if ~isempty(ALIGN.volBase.talInfo) %load talInfo if it exists
  talInfo = ALIGN.volBase.talInfo;
  talInfo.filename = ALIGN.volumePath;
  talInfo = talairach(talInfo);
elseif ~isempty(ALIGN.volBase.vol2tal) 
  % if talInfo doesn't exist, but there's a talXform defined, can use that
  % to initialize the 8 defining points (could happen if had defined a talXform
  % in another program or before the talInfo field was instituted)
  % convert vol2tal into 8 chosen points by reversing, using the pre-defined
  % talairach values of the defined points:
  points = pinv(ALIGN.volBase.vol2tal)*talPoints;  
  talInfo.AC  = points(1:3,1);
  talInfo.PC  = points(1:3,2);
  talInfo.SAC = points(1:3,3);  
  talInfo.IAC = points(1:3,4);  
  talInfo.PPC = points(1:3,5); 
  talInfo.AAC = points(1:3,6);
  talInfo.LAC = points(1:3,7);  
  talInfo.RAC = points(1:3,8);  
  clear points;
  
  talInfo.filename =  ALIGN.volumePath;
  talInfo.vol2view = eye(4);
  talInfo = talairach(talInfo);
else % if can't initialize, start from scratch
  talInfo = talairach(ALIGN.volumePath);
end

%if user hits ok, convert to a vol2tal, and save
if ~isempty(talInfo)
  points = [talInfo.AC;talInfo.PC;talInfo.SAC;talInfo.IAC;talInfo.PPC;talInfo.AAC;talInfo.LAC;talInfo.RAC]';
  points(4, :) = ones(1,size(points,2));
  talTransform = talPoints*pinv(points);
  ALIGN.volBase.vol2tal = talTransform;
  ALIGN.volBase.talInfo = talInfo;
  
  % save the matFile for the base, to save vol2tal and talInfo
  [tf base] = isbase(ALIGN.volBase);
  matFilename = sprintf('%s.mat',stripext(ALIGN.volumePath));
  base.data = [];base.hdr = []; 
  eval(sprintf('save %s base',matFilename));
  clear base 
  
  % once the vol2tal Xform has been defined, allow subjects to export it
  set(handles.exportTal2Session,'Enable','on');
end

% We've decided to keep everything in the base structure
% and not to change the NIFTI headers, e.g., not re-set
% sform_code to 3, and not change s-form to the TalXform.
% Rather, leave sform_code as 1, leave sform as alignment,
% and save talXform to the base structure.

% --------------------------------------------------------------------
function exportTal2Session_Callback(hObject, eventdata, handles)
% hObject    handle to ExportTal2Session (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global ALIGN

talInfo = ALIGN.volBase.talInfo;
vol2tal = ALIGN.volBase.vol2tal;
vol2mag = ALIGN.volBase.vol2mag;

exportTal2mrLR(vol2tal, vol2mag, talInfo);  
  


% --- Executes on button press in Rot1L.
function Rot1L_Callback(hObject, eventdata, handles)
% hObject    handle to Rot1L (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
incrementRot(handles,1,-1);

% --- Executes on button press in Rot1R.
function Rot1R_Callback(hObject, eventdata, handles)
% hObject    handle to Rot1R (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
incrementRot(handles,1,1);

% --- Executes on button press in Rot2L.
function Rot2L_Callback(hObject, eventdata, handles)
% hObject    handle to Rot2L (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
incrementRot(handles,2,-1);

% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
incrementRot(handles,2,1);

% --- Executes on button press in Rot3L.
function Rot3L_Callback(hObject, eventdata, handles)
% hObject    handle to Rot3L (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
incrementRot(handles,3,-1);

% --- Executes on button press in Rot3R.
function Rot3R_Callback(hObject, eventdata, handles)
% hObject    handle to Rot3R (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
incrementRot(handles,3,1);

% --- Executes on button press in Trans1L.
function Trans1L_Callback(hObject, eventdata, handles)
% hObject    handle to Trans1L (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
incrementTrans(handles,1,-1);

% --- Executes on button press in Trans1R.
function Trans1R_Callback(hObject, eventdata, handles)
% hObject    handle to Trans1R (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
incrementTrans(handles,1,1);

% --- Executes on button press in Trans2L.
function Trans2L_Callback(hObject, eventdata, handles)
% hObject    handle to Trans2L (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
incrementTrans(handles,2,-1);

% --- Executes on button press in Trans2R.
function Trans2R_Callback(hObject, eventdata, handles)
% hObject    handle to Trans2R (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
incrementTrans(handles,2,1);

% --- Executes on button press in Trans3L.
function Trans3L_Callback(hObject, eventdata, handles)
% hObject    handle to Trans3L (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
incrementTrans(handles,3,-1);

% --- Executes on button press in Trans3R.
function Trans3R_Callback(hObject, eventdata, handles)
% hObject    handle to Trans3R (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
incrementTrans(handles,3,1);

%%%%%%%%%%%%%%%%%%%%%%%%
%    incrementTrans   %%
%%%%%%%%%%%%%%%%%%%%%%%%
function incrementTrans(handles,axisnum,incdec)

global ALIGN
[trans rot] = getAlignGUI(handles);
% command/shift makes for a smaller/larger increment
fignum = handles.figure1;
inc = 1;
if ~isempty(fignum) 
  if any(strcmp(get(fignum,'CurrentModifier'),'shift'))
    inc = 10;
  elseif any(strcmp(get(fignum,'CurrentModifier'),'command'))
    inc = 0.1;
  end
end
trans(axisnum) = trans(axisnum)+incdec*inc;
setAlignGUI(handles,'trans',trans);
ALIGN.guiXform = getGuiXform(handles);
refreshAlignDisplay(handles);

%%%%%%%%%%%%%%%%%%%%%%
%     incrementRot  %%
%%%%%%%%%%%%%%%%%%%%%%
function incrementRot(handles,axisnum,incdec)

global ALIGN
[trans rot] = getAlignGUI(handles);
% command/shift makes for a smaller/larger increment
fignum = handles.figure1;
inc = 1;
if ~isempty(fignum)
  if any(strcmp(get(fignum,'CurrentModifier'),'shift'))
    inc = 10;
  elseif any(strcmp(get(fignum,'CurrentModifier'),'command'))
    inc = 0.1;
  end
end
rot(axisnum) = rot(axisnum)+incdec*inc;
setAlignGUI(handles,'rot',rot);
ALIGN.guiXform = getGuiXform(handles);
refreshAlignDisplay(handles);

%%%%%%%%%%%%%%%%%%%%%%%
%     setAlignTitle  %%
%%%%%%%%%%%%%%%%%%%%%%%
function setAlignTitle(handles)

global ALIGN

if ~isfield(ALIGN,'volumePath')  || isempty(ALIGN.volumePath)
  volumeName = '(none)';
else
  volumeName = getLastDir(ALIGN.volumePath);
end

if ~isfield(ALIGN,'inplanePath') || isempty(ALIGN.inplanePath)
  inplaneName = '(none)';
else
  inplaneName = getLastDir(ALIGN.inplanePath);
end

set(handles.figure1,'name',sprintf('mrAlign: %s -> %s',inplaneName,volumeName));


