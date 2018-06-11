function varargout = ea_waveletfusion(varargin)
% EA_WAVELETFUSION MATLAB code for ea_waveletfusion.fig
%      EA_WAVELETFUSION, by itself, creates a new EA_WAVELETFUSION or raises the existing
%      singleton*.
%
%      H = EA_WAVELETFUSION returns the handle to a new EA_WAVELETFUSION or the handle to
%      the existing singleton*.
%
%      EA_WAVELETFUSION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in EA_WAVELETFUSION.M with the given input arguments.
%
%      EA_WAVELETFUSION('Property','Value',...) creates a new EA_WAVELETFUSION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ea_waveletfusion_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ea_waveletfusion_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ea_waveletfusion

% Last Modified by GUIDE v2.5 31-Jul-2017 20:12:55

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ea_waveletfusion_OpeningFcn, ...
                   'gui_OutputFcn',  @ea_waveletfusion_OutputFcn, ...
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


% --- Executes just before ea_waveletfusion is made visible.
function ea_waveletfusion_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ea_waveletfusion (see VARARGIN)

% Choose default command line output for ea_waveletfusion
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

earoot=ea_getearoot;
im=imread([earoot,'icons',filesep,'logo_lead_dbs.png']);
image(im);
axis off;
axis equal;
% UIWAIT makes ea_waveletfusion wait for user response (see UIRESUME)
% uiwait(handles.ea_waveletfusion);


% --- Outputs from this function are returned to the command line.
function varargout = ea_waveletfusion_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in axial.
function axial_Callback(hObject, eventdata, handles)
% hObject    handle to axial (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[fi,pth]=uigetfile('*.nii','Choose axial volume...');
if ~fi
    set(hObject,'String','Axial Acquisition');
    setappdata(handles.ea_waveletfusion,'axial',[]);
else
    set(hObject,'String',fullfile(pth,fi));
    setappdata(handles.ea_waveletfusion,'axial',fullfile(pth,fi));
end
% --- Executes on button press in coronal.
function coronal_Callback(hObject, eventdata, handles)
% hObject    handle to coronal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[fi,pth]=uigetfile('*.nii','Choose coronal volume...');
if ~fi
    set(hObject,'String','Coeonal Acquisition');
    setappdata(handles.ea_waveletfusion,'coronal',[]);
else

set(hObject,'String',fullfile(pth,fi));
setappdata(handles.ea_waveletfusion,'coronal',fullfile(pth,fi));
end

% --- Executes on button press in sagittal.
function sagittal_Callback(hObject, eventdata, handles)
% hObject    handle to sagittal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[fi,pth]=uigetfile('*.nii','Choose sagittal volume...');
if ~fi
    set(hObject,'String','Sagittal Acquisition');
    setappdata(handles.ea_waveletfusion,'sagittal',[]);
else
    set(hObject,'String',fullfile(pth,fi));
    setappdata(handles.ea_waveletfusion,'sagittal',fullfile(pth,fi));
end

% --- Executes on button press in fusevolumes.
function fusevolumes_Callback(hObject, eventdata, handles)
% hObject    handle to fusevolumes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[fi,pth]=uiputfile('*.nii','Save fused volume as...','fused.nii');
ea_busyaction('on',handles.ea_waveletfusion,'wavelet');
axial=getappdata(handles.ea_waveletfusion,'axial');
coronal=getappdata(handles.ea_waveletfusion,'coronal');
sagittal=getappdata(handles.ea_waveletfusion,'sagittal');

if ~isempty(axial)
    axnii=ea_load_nii(axial);
    axv=axnii.img;
else
    axv=[];
end
if ~isempty(coronal)
    cornii=ea_load_nii(coronal);
    corv=cornii.img;
else
    corv=[];
end
if ~isempty(sagittal)
    sagnii=ea_load_nii(sagittal);
    sagv=sagnii.img;
else
    sagv=[];
end

fusev=fuseVolW(axv,corv,sagv);

fusenii=axnii;
fusenii.img=fusev;
fusenii.fname=fullfile(pth,fi);
ea_write_nii(fusenii);
ea_busyaction('off',handles.ea_waveletfusion,'wavelet');
ea_chirp
ea_methods('',...
            ['Multiple acquisitions were fused to a combined volume based on a wavelet-based 3D image fusion approach ',...
            '(Aganj 2012).'],...
            {'I. Aganj, C. Lenglet, E. Yacoub, G. Sapiro, and N. Harel, ?A 3D wavelet fusion approach for the reconstruction of isotropic-resolution MR images from orthogonal anisotropic-resolution scans,? Magnetic Resonance in Medicine, vol. 67, no. 4, pp. 1167?1172, 2012.'});
msgbox('Successfully exported fused volume.','Success');
