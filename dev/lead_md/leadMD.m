function varargout = leadMD(varargin)
% LEADMD MATLAB code for leadMD.fig
%      LEADMD, by itself, creates a new LEADMD or raises the existing
%      singleton*.
%
%      H = LEADMD returns the handle to a new LEADMD or the handle to
%      the existing singleton*.
%
%      LEADMD('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in LEADMD.M with the given input arguments.
%
%      LEADMD('Property','Value',...) creates a new LEADMD or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before leadMD_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to leadMD_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help leadMD

% Last Modified by GUIDE v2.5 31-Aug-2015 22:07:26

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @leadMD_OpeningFcn, ...
    'gui_OutputFcn',  @leadMD_OutputFcn, ...
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


% --- Executes just before leadMD is made visible.
function leadMD_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to leadMD (see VARARGIN)

% Choose default command line output for leadMD
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);


set(handles.versiontxt,'String',['v',ea_getvsn('local')]);

% get electrode model specs and place in popup
set(handles.electrode_model_popup,'String',ea_resolve_elspec);

set(gcf,'name','Welcome to LEAD-DBS MD');

im=imread('ea_logo.png');
image(im);
axis off;
axis equal;

% UIWAIT makes leadMD wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = leadMD_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in impdatatabbut.
function impdatatabbut_Callback(hObject, eventdata, handles)
% hObject    handle to impdatatabbut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in reviewtabbut.
function reviewtabbut_Callback(hObject, eventdata, handles)
% hObject    handle to reviewtabbut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in viztabbut.
function viztabbut_Callback(hObject, eventdata, handles)
% hObject    handle to viztabbut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in dtiax.
function dtiax_Callback(hObject, eventdata, handles)
% hObject    handle to dtiax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
options.prefs=ea_prefs('tmp');
[img,path]=ea_importdata(options.prefs.dti);
hObject.CData=img;
hObject.String='';
setappdata(hObject,'path',path);

% --- Executes on button press in preopax.
function preopax_Callback(hObject, eventdata, handles)
% hObject    handle to preopax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
options.prefs=ea_prefs('tmp');
[img,path]=ea_importdata(options.prefs.prenii_unnormalized);
hObject.CData=img;
hObject.String='';
setappdata(hObject,'path',path);

% --- Executes on button press in postopax.
function postopax_Callback(hObject, eventdata, handles)
% hObject    handle to postopax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
options.prefs=ea_prefs('tmp');
[img,path]=ea_importdata(options.prefs.tranii_unnormalized);
hObject.CData=img;
hObject.String='';
setappdata(hObject,'path',path);


% --- Executes on button press in isalreadycoregisted.
function isalreadycoregisted_Callback(hObject, eventdata, handles)
% hObject    handle to isalreadycoregisted (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of isalreadycoregisted


% --- Executes on selection change in MRCT.
function MRCT_Callback(hObject, eventdata, handles)
% hObject    handle to MRCT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns MRCT contents as cell array
%        contents{get(hObject,'Value')} returns selected item from MRCT


% --- Executes during object creation, after setting all properties.
function MRCT_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MRCT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in corax.
function corax_Callback(hObject, eventdata, handles)
% hObject    handle to corax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
options.prefs=ea_prefs('tmp');
[img,path]=ea_importdata(options.prefs.cornii_unnormalized);
hObject.CData=img;
hObject.String='';
setappdata(hObject,'path',path);

% --- Executes on button press in sagax.
function sagax_Callback(hObject, eventdata, handles)
% hObject    handle to sagax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
options.prefs=ea_prefs('tmp');
[img,path]=ea_importdata(options.prefs.sagnii_unnormalized);
hObject.CData=img;
hObject.String='';
setappdata(hObject,'path',path);

% --- Executes on button press in runbutton.
function runbutton_Callback(hObject, eventdata, handles)
% hObject    handle to runbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

keyboard

outfolder=uigetdir('','Specify output directory...');
if outfolder

    %% run trajectory reconstruction.
    options=ea_step1options(handles);
    options.root=[fileparts(fileparts(outfolder)),filesep];
    [~,options.patientname]=fileparts(fileparts(outfolder));
    keyboard

    % copy files from temp to new folder
    fis={options.prefs.tranii};

    ea_autocoord(options);

end






function [img,opath]=ea_importdata(fname)

fspec{1,1}='*.nii'; fspec{2,1}='*';
fspec{1,2}='Nifti-1'; fspec{2,2}='DICOM';

[fis,path,fselect]=uigetfile(fspec,'Choose .nii or DICOM data...','MultiSelect','on');

directory=[tempdir,'leadMD'];
try
    rmdir(directory,'s');
end
mkdir(directory);

options.prefs=ea_prefs('temp');
if fselect==1 % nifti
    if iscell(fis)
        fis=fis{1};
    end
    copyfile([path,fis],[directory,'import.nii']);
else % DICOM

    fips=cellfun(@(x) strcat(path,x), fis, 'UniformOutput', 0);

    matlabbatch{1}.spm.util.import.dicom.data = fips';
    matlabbatch{1}.spm.util.import.dicom.root = 'flat';
    matlabbatch{1}.spm.util.import.dicom.outdir = directory;
    matlabbatch{1}.spm.util.import.dicom.protfilter = '.*';
    matlabbatch{1}.spm.util.import.dicom.convopts.format = 'nii';
    matlabbatch{1}.spm.util.import.dicom.convopts.icedims = 0;
    cfg_util('run',{matlabbatch});
    clear matlabbatch

    di=dir([directory,'*.nii']);
    movefile([directory,di(1).name],[directory,'import.nii']);
end

opath=[directory,options.prefs.dti];
movefile([directory,'import.nii'],opath);
nii=nifti(opath);
img=squeeze(nii.dat(round(size(nii.dat,1)/2),:,:));


img=(img-min(img(:)))/(max(img(:))-min(img(:)));
img=repmat(img,1,1,3);


function options=ea_step1options(handles)


%% some manual options that can be set:


options.endtolerance=10; % how many slices to use with zero signal until end of electrode estimate.
options.sprungwert=4; % how far electrode centroid may be (in xy axis) from last to current slice.
options.refinesteps=0; % how often to re-iterate to reconstruct trajectory. More than 2 should usually not be beneficial. Use 0 to use the direct measurement.
options.tra_stdfactor=0.9; % Default: 0.9 - the lower this factor, the lower the threshold (more included pixels in tra process).
options.cor_stdfactor=1.0; % Default: 1.0 - the higher this factor, the lower the threshold (more included pixels in cor process).




%% set options

%uipatdir=get(handles.patdir_choosebox,'String');

options.earoot=[fileparts(which('lead')),filesep];
options.dicomimp=0;

options.normalize.do=1;
options.normalize.method='SPM12 DARTEL nonlinear [MR/CT]';
options.normalize.methodn=1;

options.normalize.check=0;


% set modality (MR/CT) in options
options.modality = get(handles.MRCT,'Value');

if options.modality==2; % CT

    % coreg CT
    options.coregct.do=~get(handles.isalreadycoregisted,'Value');
    options.coregct.method='Coregister postop-CT with preop-MRI (Matlab Imreg)';
    options.coregct.methodn=1;
    options.coregct.coregthreshs=[0.6,0.4];

    options.coregctcheck=0;
else
    options.coregct.do=0;
end







options.verbose=3; % 4: Show figures but close them 3: Show all but close all figs except resultfig 2: Show all and leave figs open, 1: Show displays only, 0: Show no feedback.
options.sides=1:2;

options.doreconstruction=1;
options.maskwindow=10; % initialize at 10
options.automask=1; % set automask flag
options.autoimprove=0; % if true, there will be some pauses at critical points so that the process can be better visualized. Mainly for demonstration or debugging problems.

options.axiscontrast=9; % if 8: use tra only but smooth it before. % if 9: use mean of cor and tra but smooth it. % if 10: use raw tra only.
options.zresolution=10; % voxels are being parcellated into this amount of portions.

options.atl.genpt=0; % generate patient specific atlases
options.atl.normalize=0; % normalize patient specific atlasset.
options.atl.can=1; % display canonical atlases
options.atl.pt=0; % display patient specific atlases

options.manualheightcorrection=0;


options.d2.write=0;
options.d2.atlasopacity=0.15;


options.d3.write=0;
options.d3.prolong_electrode=2;
options.d3.verbose='on';
options.d3.elrendering=1;
options.d3.hlactivecontacts=0;
options.d3.showactivecontacts=1;
options.d3.showpassivecontacts=1;
options.d3.showisovolume=0;
options.d3.isovscloud=0;
options.d3.autoserver=0;

options.numcontacts=4;
options.entrypoint='STN, GPi or ViM';
options.entrypointn=1;

options.writeoutpm=0;

options.elmodeln = get(handles.electrode_model_popup,'Value');
string_list = get(handles.electrode_model_popup,'String');
options.elmodel=string_list{options.elmodeln};
options.atlasset=get(handles.atlassetpopup,'String'); %{get(handles.atlassetpopup,'Value')}
options.atlasset=options.atlasset{get(handles.atlassetpopup,'Value')};
options.atlassetn=get(handles.atlassetpopup,'Value');

if strcmp(options.atlasset,'Use none');
    options.d3.writeatlases=0;
    options.d2.writeatlases=1;
else
    options.d3.writeatlases=1;
    options.d2.writeatlases=1;
end


options.expstatvat.do=0;

options.fiberthresh=10;
options.writeoutstats=1;


options.colormap=colormap;

options.dolc=0;


% --- Executes on selection change in electrode_model_popup.
function electrode_model_popup_Callback(hObject, eventdata, handles)
% hObject    handle to electrode_model_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns electrode_model_popup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from electrode_model_popup


% --- Executes during object creation, after setting all properties.
function electrode_model_popup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to electrode_model_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
