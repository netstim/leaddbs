function varargout = uw_overlay_bwgui(varargin)
% UW_OVERLAY_BWGUI M-file for uw_overlay_bwgui.fig
% Version 1.1
%
%   GUI - uw_overlay_bwgui
%
%   An interface desiged to create and preview overlays using the
%   uw_overlay_bwthresh function.  All of the settings avaliable on the
%   command line of bwthres are avaliable inside the GUI.
%
%   You have the option to save the NIfTI file created by bwthresh, or
%   convert it to a seris of DICOM images for upload to a PACS or
%   neuroimaging system.  Specify a DICOM file from the original anatomical
%   scan, and specify a new series number and description, and
%   convert2dicom will create a new directory with DICOM files that contain
%   the original patient and scanner information.
%
%   For help on the underlying functions used by this GUI, type:
%   doc uw_overlay_bwthresh
%   doc uw_overlay_convert2dicom
%
%   Copyright 2007 Samuel A. Hurley - samuel.hurley[at]gmail.com
%   Universtiy of Wisconsin - Madison | Applied Neuro fMRI Lab
%   Published under GNU General Public License Version 2 or 3, see <http://www.gnu.org/licenses/>


%     This program is free software; you can redistribute it and/or modify
%     it under the terms of the GNU General Public License version 2 or 3
%     as published by the Free Software Foundation
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.


% Last Modified by GUIDE v2.5 19-Aug-2007 17:17:24

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @uw_overlay_bwgui_OpeningFcn, ...
                   'gui_OutputFcn',  @uw_overlay_bwgui_OutputFcn, ...
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


% --- Executes just before uw_overlay_bwgui is made visible.
function uw_overlay_bwgui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to uw_overlay_bwgui (see VARARGIN)

% Choose default command line output for uw_overlay_bwgui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes uw_overlay_bwgui wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% Define global variables used by uw_overlay_bwthresh
% Program Variables
global A; %  Handle to original anatomical image
global C;  % Anatomy image matrix
global tempOutput;  % Temporary output image handle
global overlayScaled; % Overlay image matrix, scaled to anatomy size
global anatMaxval;
global anatMinval;

% User Input variables
global anat_img; %#ok<NUSED>
global overlay_img; %#ok<NUSED>
global output_img;
global overlayThresh;
global overlayPercent;
global overlayAlpha;
global overlayInterp;
global overlayCoreg;

global coregImage;
global includeBrainVolume;

% Set default values
output_img = 'combined.img';
overlayThresh = 0;
overlayPercent = 100;
overlayAlpha = 100;
overlayInterp = 1;
overlayCoreg = eye(4);
coregImage = '0';
includeBrainVolume = false;

% Define a global switch to ensure that files are specified before the
% controls are enabled
global anatSpecified;
global overlaySpecified;
global dicomSpecified;
anatSpecified = 0;
overlaySpecified = 0;
dicomSpecified = 0;

% Print out a friendly message to the user
disp('UW Overlay B&W GUI, Version 1.0');
disp('Copyright 2007 Samuel A. Hurley - samuel.hurley[at]gmail.com');
disp('Universtiy of Wisconsin - Madison | Applied Neuro fMRI Lab');
disp('Published under GNU General Public License Version 2 or 3, see <http://www.gnu.org/licenses/>');

% Create an SPM graphics window
	fg = spm_figure('Findwin','Graphics');
	if isempty(fg),
		fg=spm_figure('Create','Graphics');
		if isempty(fg),
			error('Cant create graphics window');
		end;
    else
		spm_figure('Clear','Graphics');
	end;

% --- Outputs from this function are returned to the command line.
function varargout = uw_overlay_bwgui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in radiobutton1.
function radiobutton1_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton1

global overlayCoreg;
global anat_img;
global overlay_img;

overlayCoreg = eye(4); % Reset coregistration matrix back to original
set(handles.edit2, 'String', 'not done');

set(handles.radiobutton1, 'Value', 1);
set(handles.radiobutton2, 'Value', 0);

set(handles.edit1, 'Enable', 'off');
set(handles.edit2, 'Enable', 'off');
set(handles.pushbutton1, 'Enable', 'off');
set(handles.pushbutton2, 'Enable', 'off');

rescaleOverlay(anat_img, overlay_img, handles);  % Need to re-slice using original eye(4) transformation


% --- Executes on button press in radiobutton2.
function radiobutton2_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton2

global coregImage;

set(handles.radiobutton1, 'Value', 0);
set(handles.radiobutton2, 'Value', 1);

set(handles.edit1, 'Enable', 'inactive');
set(handles.edit2, 'Enable', 'inactive');
set(handles.pushbutton1, 'Enable', 'on');

if coregImage ~= '0'
    set(handles.pushbutton2, 'Enable', 'on');   % Only enable the calculate button if
                                                % an input image has been specified
end

function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global coregImage;

[filename status] = spm_select(1, 'image', 'Select source image...');
if status == 1
    coregImage = filename;
    set(handles.edit1, 'String', filename);
    set(handles.pushbutton2, 'Enable', 'on');
end


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global anat_img;
global overlay_img;

global coregImage;
global overlayCoreg;


k = waitbar(.2, 'Loading images...');
% Load the two images
A = spm_vol(anat_img);
waitbar(.3, k);
B = spm_vol(coregImage);

waitbar(.4, k, 'Coregistering. Go get a cup of tea, this may take a while...');
% Figure out the transform
m = spm_coreg(A, B);

waitbar(.6, k, 'Reslicing...');
% Convert using spm_matrix
overlayCoreg = spm_matrix(m);


% Re-run the rescale operation on the two images, since the value of the
% transform matrix M has now changed
waitbar(.7, k);
rescaleOverlay(anat_img, overlay_img, handles);

% Tell the user we're FINALLY done.
delete(k);
set(handles.edit2, 'String', 'Done.');

function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on text3 modification
function text3_Callback(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

global overlayThresh;
overlayThresh = str2num(get(hObject, 'String'));

if overlayThresh > -50 && overlayThresh < 150
    set(handles.slider1, 'Value', overlayThresh)
end

setPreviewButton(handles);

% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

global overlayThresh;

overlayThresh = round(get(hObject, 'Value'));   % Only accept integer values
set(handles.text3, 'String', num2str(overlayThresh));
setPreviewButton(handles);


% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider2_Callback(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

global overlayPercent;

overlayPercent = round(get(hObject, 'Value'));
set(handles.text5, 'String', [num2str(overlayPercent) '%']);
setPreviewButton(handles);

% --- Executes during object creation, after setting all properties.
function slider2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider3_Callback(hObject, eventdata, handles)
% hObject    handle to slider3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

global overlayAlpha;

overlayAlpha = round(get(hObject, 'Value'));
set(handles.text7, 'String', [num2str(overlayAlpha) '%']);
setPreviewButton(handles);

% --- Executes during object creation, after setting all properties.
function slider3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% Uncomment to see whats going on!
% global anat_img;
% global overlay_img;
% global overlayThresh;
% global overlayPercent;
% global overlayAlpha;
% global overlayInterp;
% global overlayCoreg;

% anat_img
% overlay_img
% overlayThresh
% overlayPercent
% overlayAlpha
% overlayInterp
% overlayCoreg

% previewOverlay is a faster way to do this:
% uw_overlay_bwthresh(anat_img, overlay_img, 'temp.img', overlayThresh, overlayPercent, overlayAlpha, overlayInterp, overlayCoreg);
% spm_check_registration temp.img

previewOverlay();
resetPreviewButton(handles);

% --- Executes on button press in radiobutton3.
function radiobutton3_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton3

set(handles.radiobutton4, 'Value', 0);
set(handles.radiobutton3, 'Value', 1);

set(handles.edit3, 'Enable', 'on');
set(handles.pushbutton4, 'Enable', 'on');

set(handles.pushbutton7, 'Enable', 'off');
set(handles.pushbutton8, 'Enable', 'off');
set(handles.edit8, 'Enable', 'off');
set(handles.edit9, 'Enable', 'off');


% --- Executes on button press in radiobutton4.
function radiobutton4_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton4

global dicomSpecified;

set(handles.radiobutton4, 'Value', 1);
set(handles.radiobutton3, 'Value', 0);

set(handles.edit3, 'Enable', 'off');
set(handles.pushbutton4, 'Enable', 'off');
if dicomSpecified ==1
    set(handles.pushbutton7, 'Enable', 'on');
end
set(handles.pushbutton8, 'Enable', 'on');
set(handles.edit8, 'Enable', 'on');
set(handles.edit9, 'Enable', 'on');

function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double


% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

k = waitbar(.2, 'Saving image...');
newFileName = deblank(get(handles.edit3, 'String'));
copyfile('temp.img', [newFileName, '.img']);
k = waitbar(.4, k, 'Saving header...');
copyfile('temp.hdr', [newFileName, '.hdr']);
waitbar(.6, k, 'Removing temp files...');
delete('temp.img');
delete('temp.hdr');
waitbar(.8, k, 'Displaying file...');

spm_check_registration([newFileName, '.img']);
if exist('k'); delete(k); end;

function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double


% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit5_Callback(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit5 as text
%        str2double(get(hObject,'String')) returns contents of edit5 as a double


% --- Executes during object creation, after setting all properties.
function edit5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global anat_img;
global overlay_img;
global anatSpecified;
global overlaySpecified;

[filename status] = spm_select(1, 'image', 'Select anatomical...');
if status == 1
    anat_img = filename;
    anatSpecified = 1;
    set(handles.edit4, 'String', filename);
end

if anatSpecified == 1 & overlaySpecified == 1
    set(handles.radiobutton1, 'Enable', 'on');
    set(handles.radiobutton2, 'Enable', 'on');
    %set(handles.radiobutton3, 'Enable', 'on');  % Make program more idiot-proof
    %set(handles.radiobutton4, 'Enable', 'on');
    set(handles.slider1, 'Enable', 'on');
    set(handles.slider2, 'Enable', 'on');
    set(handles.slider3, 'Enable', 'on');
    set(handles.txtInterp, 'Enable', 'on');
    set(handles.pushbutton3, 'Enable', 'on');
    %set(handles.pushbutton4, 'Enable', 'on');
    set(handles.text3, 'Enable', 'on');
    set(handles.edit3, 'Enable', 'on');
    
    rescaleOverlay(anat_img, overlay_img, handles);  % Function to load the two images and
                                            % then rescale the overlay to
                                            % the same size as anatomy
end
    

% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global anat_img;
global overlay_img;
global anatSpecified;
global overlaySpecified;

global overlayInterp;

[filename status] = spm_select(1, 'image', 'Select overlay...');
if status == 1
    overlay_img = filename;
    overlaySpecified = 1;
    set(handles.edit5, 'String', filename);
end

if anatSpecified == 1 & overlaySpecified == 1
    set(handles.radiobutton1, 'Enable', 'on');
    set(handles.radiobutton2, 'Enable', 'on');
    %set(handles.radiobutton3, 'Enable', 'on');     % Make program more idiot-proof
    %set(handles.radiobutton4, 'Enable', 'on');
    set(handles.slider1, 'Enable', 'on');
    set(handles.slider2, 'Enable', 'on');
    set(handles.slider3, 'Enable', 'on');
    set(handles.txtInterp, 'Enable', 'on');
    set(handles.pushbutton3, 'Enable', 'on');
    %set(handles.pushbutton4, 'Enable', 'on');
    set(handles.text3, 'Enable', 'on');
    set(handles.edit3, 'Enable', 'on');

    rescaleOverlay(anat_img, overlay_img, handles);  % Function to load the two images and
                                            % then rescale the overlay to
                                            % the same size as anatomy
end

function txtInterp_Callback(hObject, eventdata, handles)
% hObject    handle to txtInterp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtInterp as text
%        str2double(get(hObject,'String')) returns contents of txtInterp as a double

global overlayInterp;
global anat_img;
global overlay_img;

overlayInterp = round(str2num(get(hObject, 'String')));

rescaleOverlay(anat_img, overlay_img, handles);

% --- Executes during object creation, after setting all properties.
function txtInterp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtInterp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% <-- Function to load the anatomy and overlay into A and B, then rescale B
function rescaleOverlay(anat_img, overlay_img, handles)

global A;
global C;
global overlayScaled;
global anatMaxval;
global anatMinval;
global overlayInterp;

global overlayCoreg; % If the user later specifies a coregistration image,
                     % This function rescaleOverlay will need to be re-run
                     % with the new overlayCoreg

k = waitbar(0, 'Loading images...');
% Load images with spm_vol
A = spm_vol(anat_img);      % A - anatomy image datastructure
B = spm_vol(overlay_img);   % B - overlay image datastructure 

waitbar(.1, k, 'Reading images...');
C = spm_read_vols(A);       % C - anatomy image matrix

anatDim = A.dim;            % dimensions of anatomy image matrix
anatMat = A.mat;            % anatomy transformation matrix
overlayMat = B.mat;         % overlay transformation matrix

% Initialize other variables ahead of time
M = zeros(4);               % M - transformation matrix to scale up overlay onto anatomy
overlayScaled = zeros(anatDim);  % A scaled up version of overlay data

% Find maximum and minimum values in anat. image matrix C
waitbar(.2, k, 'Finding max value...');
elementsInC = anatDim(1)*anatDim(2)*anatDim(3);
Creshape    = reshape(C, 1, elementsInC);       % Use reshape to convert C into a vector

anatMaxval  = max(Creshape);  % anatMaxval - maximum value in anatomy image
anatMinval  = min(Creshape);  % anatMinVal - minimum value in anatomy image


% IV: Resample the overlay B to the same dimesions as the anat image C, using transform matrix M
voxelM = overlayMat \ overlayCoreg * anatMat;     % Transforms from voxel coords in anatomy to voxel coords in overlay
sliceM = 0;

waitbar(.3, k, 'Resizing overlay...')
for n = 1:anatDim(3)   % Loop through and reslice overlay image
    sliceM(3) = n;
    M = voxelM * spm_matrix(sliceM);  % First shift to the nth slice in anat coords, then transform to overlay coords
    overlayScaled(:,:,n) = spm_slice_vol(B, M, anatDim(1:2), overlayInterp);
    waitbar(.3 + (n/anatDim(3) * .7), k);
end
delete(k);
setPreviewButton(handles);

% <-- Function to calculate and plot a preview image in the SPM graphics window
function previewOverlay()

global A;
global C;
global overlayScaled;
global anatMaxval;
global anatMinval;

global overlayThresh;
global overlayPercent;
global overlayAlpha;

global includeBrainVolume;

global anat_img;
global tempOutput;

k = waitbar(0, 'Computing overlay brightness...');
% Determine the 'brightness' of the overlay, with specified overlayPercent
overlayPixelvalue = (overlayPercent/100) * (anatMaxval - anatMinval) + anatMinval;

waitbar(.1, k, 'Computing overlay mask...');
% V:  Calculate masks for anatomy and overlay image
overlayMask = overlayScaled >= overlayThresh; % Mask where overlay is greater than threshold
                                              % 0 for regions below thresh, 1 for regions at or
                                              % above threshold
                                              

% If the overlay lies outside the brain/anatomy image volume, DON'T overlay
% it!!
% In other words, only place an overlay into  non-zero voxels of the
% original images:
if includeBrainVolume == true
    overlayMask = overlayMask .* (C ~= 0);
end

anatMask = C .* (~overlayMask);  % Use overlay mask to remove elements from anat. image

waitbar(.2, k, 'Combining images...');
% VI: Combine overlayMask and anatMask with overlayPixelvalue and overlayAlpha parameters
alphaBlend = round((C .* (1 - (overlayAlpha/100)) + overlayMask * overlayPixelvalue * (overlayAlpha/100)));

waitbar(.3, k);
D = anatMask + alphaBlend .* overlayMask;       % D is now our finished image matrix
D = int16(D);


waitbar(.4, k, 'Writing new image...');
% VII: Write to a temporary output file
A.fname = 'temp.img';   % Then change the name to a temp. file
A.descrip = 'Temporary output from uw_overlay_bwgui.m';

% Write out the temp file
% Don't use spm_write_vol, because it rescales the data values:
% A = spm_write_vol(A, D);    % Write the temp file

waitbar(.45, k);
A = spm_create_vol(A);  % Create new header
waitbar(.5, k);

n = size(D, 3);
for m = 1:n
    waitbar(m/n*.4 + .5, k);
    A = spm_write_plane(A, D(:,:,m), m);  % Create new image
end

waitbar(.9, k, 'Displaying preview...');
spm_check_registration temp.img;
delete(k);


function edit7_Callback(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit7 as text
%        str2double(get(hObject,'String')) returns contents of edit7 as a double


% --- Executes during object creation, after setting all properties.
function edit7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit8_Callback(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit8 as text
%        str2double(get(hObject,'String')) returns contents of edit8 as a double


% --- Executes during object creation, after setting all properties.
function edit8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit9_Callback(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit9 as text
%        str2double(get(hObject,'String')) returns contents of edit9 as a double


% --- Executes during object creation, after setting all properties.
function edit9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global dicomSpecified;

if dicomSpecified == 1   % Check that the DICOM file has been specified
    dicom_file = get(handles.edit7, 'String');
    merged_file = 'temp.img';
    newSeriesNumber = str2num(get(handles.edit8, 'String'));
    newSeriesDescription = get(handles.edit9, 'String');
    outputDirectory = num2str(newSeriesNumber);
    mergedImageVolume = 1;
    outputImagePosition = 2;
    uw_overlay_convert2dicom(dicom_file, merged_file, newSeriesNumber, newSeriesDescription, outputDirectory, mergedImageVolume, outputImagePosition);
    waitbar(1, 'DICOM output complete!');
else
    waitbar(1, 'Error, DICOM file not specified');
end


% --- Executes on button press in pushbutton8.
function pushbutton8_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global dicomSpecified;

[dicomfile status] = spm_select(1, 'any', 'Select DICOM file...');
if status == 1
    dicomSpecified = 1;
    set(handles.pushbutton7, 'Enable', 'on');
    set(handles.edit7, 'String', dicomfile);
end

%----------
function resetPreviewButton(handles)

set(handles.pushbutton3, 'BackgroundColor', [.2 1 .2]);
set(handles.text18,      'String', ' ');

set(handles.radiobutton3, 'Enable', 'on');
set(handles.radiobutton4, 'Enable', 'on');

if get(handles.radiobutton3, 'Value') == 1
    radiobutton3_Callback(0, 0, handles)
else
    radiobutton4_Callback(0, 0, handles)
end

%-----------
function setPreviewButton(handles)

set(handles.pushbutton3, 'BackgroundColor', [1 .2 .2]);
set(handles.text18,      'String', 'Not updated ->');

set(handles.radiobutton3, 'Enable', 'off');
set(handles.radiobutton4, 'Enable', 'off');

set(handles.pushbutton4,  'Enable', 'off');
set(handles.pushbutton7,  'Enable', 'off');
set(handles.pushbutton8,  'Enable', 'off');

set(handles.edit8,        'Enable', 'off');
set(handles.edit9,        'Enable', 'off');



% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox1
global includeBrainVolume;

if get(hObject, 'Value') == 1
    includeBrainVolume = true;
else
    includeBrainVolume = false;
end