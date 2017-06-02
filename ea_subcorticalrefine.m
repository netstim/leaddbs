function varargout = ea_subcorticalrefine(varargin)
% EA_SUBCORTICALREFINE MATLAB code for ea_subcorticalrefine.fig
%      EA_SUBCORTICALREFINE, by itself, creates a new EA_SUBCORTICALREFINE or raises the existing
%      singleton*.
%
%      H = EA_SUBCORTICALREFINE returns the handle to a new EA_SUBCORTICALREFINE or the handle to
%      the existing singleton*.
%
%      EA_SUBCORTICALREFINE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in EA_SUBCORTICALREFINE.M with the given input arguments.
%
%      EA_SUBCORTICALREFINE('Property','Value',...) creates a new EA_SUBCORTICALREFINE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ea_subcorticalrefine_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ea_subcorticalrefine_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ea_subcorticalrefine

% Last Modified by GUIDE v2.5 05-Mar-2017 12:51:37

% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ea_subcorticalrefine_OpeningFcn, ...
                   'gui_OutputFcn',  @ea_subcorticalrefine_OutputFcn, ...
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


% --- Executes just before ea_subcorticalrefine is made visible.
function ea_subcorticalrefine_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ea_subcorticalrefine (see VARARGIN)




options=varargin{1};
directory=[options.root,options.patientname,filesep];

setappdata(handles.scrf,'options',options);
setappdata(handles.scrf,'directory',directory);

[pth,patientname]=fileparts(fileparts(directory));
handles.patientname.String=patientname;

set(handles.scrf,'name',['Brainshift Correction: ',patientname]);
ea_refreshscrf(options,handles,directory);
switch options.prefs.scrf.auto
    case 'nomask'
        handles.mask0.Value=1;
        handles.mask1.Value=0;
        handles.mask2.Value=0;
        %if ~exist([directory,'scrf',filesep,'scrf_instore.mat'],'file') && ~exist([directory,'scrf',filesep,'scrf.mat'],'file')
            ea_compute_scrf(handles)
        %end
    case 'mask1'
        handles.mask0.Value=0;
        handles.mask1.Value=1;
        handles.mask2.Value=0;
        %if ~exist([directory,'scrf',filesep,'scrf_instore.mat'],'file') && ~exist([directory,'scrf',filesep,'scrf.mat'],'file')
            ea_compute_scrf(handles)
        %end
    case 'mask2'
        handles.mask0.Value=0;
        handles.mask1.Value=0;
        handles.mask2.Value=1;
        %if ~exist([directory,'scrf',filesep,'scrf_instore.mat'],'file') && ~exist([directory,'scrf',filesep,'scrf.mat'],'file')
            ea_compute_scrf(handles)
        %end
end


% Choose default command line output for ea_subcorticalrefine
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ea_subcorticalrefine wait for user response (see UIRESUME)

if options.d2.write || options.d3.write
uiwait(handles.scrf);
end

function ea_refreshscrf(options,handles,directory)

standardslice=ea_loadrefineslice(directory,options,0);
refineslice=ea_loadrefineslice(directory,options,1);
set(0,'CurrentFigure',handles.scrf);


handles.scrf.CurrentAxes=handles.standardax;
imshow(standardslice);
handles.scrf.CurrentAxes=handles.scfax;
imshow(refineslice);

% display matrix:
if exist([directory,'scrf',filesep,'scrf_instore.mat'],'file')
load([directory,'scrf',filesep,'scrf_instore.mat']);
mat=ea_antsmat2mat(AffineTransform_float_3_3,fixed);
handles.affmatrix.String=sprintf('%0.2f %0.2f %0.2f %0.2f\n%0.2f %0.2f %0.2f %0.2f\n%0.2f %0.2f %0.2f %0.2f\n%0.2f %0.2f %0.2f %0.2f',mat');
end

function slice=ea_loadrefineslice(directory,options,refine)

switch refine
    case 1
        refstr='scrf';
    case 0
        refstr='standard';
end

ea_createrefineslice(directory,options,refine);
try
slice=imread([directory,'scrf',filesep,refstr,'.png']);
catch
    slice=imread([ea_getearoot,'helpers',filesep,'gui',filesep,'scrf_msg.png']);
end



function fn=ea_stripex(fn)
[~,fn]=fileparts(fn);

function ea_createrefineslice(directory,options,refine)


switch refine
    case 1
        scrf='scrf';
    case 0
        scrf='';
end


ea_createbbfiles(directory); % needs to unfortunately be done each time since coregistration may have changed.
ea_createmovim(directory,options);
ea_gencoregcheckfigs_scrf(directory,scrf,options);


function ea_createbbfiles(directory)
[options.root,options.patientname]=fileparts(fileparts(directory));
options.root=[options.root,filesep];
options.earoot=ea_getearoot;
options.prefs=ea_prefs(options.patientname);
options=ea_assignpretra(options);
if ~exist([directory,'scrf',filesep,options.prefs.prenii_unnormalized],'file')
    if ~exist([directory,'scrf'],'dir')
        mkdir([directory,'scrf']);
    end
    to{1}=[directory,'scrf',filesep,'bb.nii'];
    from{1}=[ea_space,'bb.nii'];

    ea_apply_normalization_tofile(options,from,to,[options.root,options.patientname,filesep],1);
    ea_crop_nii([directory,'scrf',filesep,'bb.nii']);

    % do put in primary anat file ? needs to be done only once.
    fis={options.prefs.prenii_unnormalized};
    copyfile([directory,fis{1}],[directory,'scrf',filesep,fis{1}])
    ea_conformspaceto([directory,'scrf',filesep,'bb.nii'],[directory,'scrf',filesep,fis{1}]);
    % cleanup:
    delete([directory,'scrf',filesep,'bb.nii']);
end
% apply tonemapping if needed
if ~exist([directory,'tp_',options.prefs.ctnii_coregistered],'file') && exist([directory,options.prefs.ctnii_coregistered],'file')
   ea_tonemapct_file(options,'native');
end
fis={options.prefs.tranii_unnormalized,options.prefs.cornii_unnormalized,options.prefs.sagnii_unnormalized,['tp_',options.prefs.ctnii_coregistered]};
for fi=1:length(fis)
    if exist([directory,fis{fi}],'file')
        copyfile([directory,fis{fi}],[directory,'scrf',filesep,fis{fi}])
        ea_conformspaceto([directory,'scrf',filesep,options.prefs.prenii_unnormalized],[directory,'scrf',filesep,fis{fi}]);
    end
end

% --- Outputs from this function are returned to the command line.
function varargout = ea_subcorticalrefine_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
%varargout{1} = handles.output;


% --- Executes on selection change in methodbutton.
function methodbutton_Callback(hObject, eventdata, handles)
% hObject    handle to methodbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns methodbutton contents as cell array
%        contents{get(hObject,'Value')} returns selected item from methodbutton


% --- Executes during object creation, after setting all properties.
function methodbutton_CreateFcn(hObject, eventdata, handles)
% hObject    handle to methodbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in computebutn.
function computebutn_Callback(hObject, eventdata, handles)
% hObject    handle to computebutn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

ea_compute_scrf(handles)

function ea_compute_scrf(handles)

options=getappdata(handles.scrf,'options');
directory=getappdata(handles.scrf,'directory');

options.coregmr.method='ANTs';

if ~handles.mask0.Value
    if ~exist([directory,'scrf',filesep,'bgmsk.nii'],'file')
        ea_addtsmask(options,1);
        for msk=1:2
            if msk==2
                btts='2';
            else
                btts='';
            end
            if ~exist([directory,'scrf',filesep],'dir')
                mkdir([directory,'scrf',filesep]);
            end
            movefile([directory,'bgmsk',btts,'.nii'],[directory,'scrf',filesep,'bgmsk',btts,'.nii']);
            ea_conformspaceto([directory,'scrf',filesep,options.prefs.prenii_unnormalized],[directory,'scrf',filesep,'bgmsk',btts,'.nii'],0);
        end
    end
    for msk=1:2
        if msk==2
            btts='2';
        else
            btts='';
        end
        msks{msk}=[directory,'scrf',filesep,'bgmsk',btts,'.nii'];
    end
end
otherfiles={};

if handles.mask1.Value
   msks=msks(1); % only use first mask.
end
if ~exist('msks','var')
    msks={};
end

ea_coreg2images(options,[directory,'scrf',filesep,'movim.nii'],[directory,'scrf',filesep,options.prefs.prenii_unnormalized],...
    [directory,'scrf',filesep,'scrfmovim.nii'],otherfiles,1,msks);

%movefile([directory,'scrf',filesep,'movim.nii'],[directory,'scrf',filesep,'scrfmovim.nii']);
movefile([directory,'scrf',filesep,'raw_movim.nii'],[directory,'scrf',filesep,'movim.nii']);

movefile([directory,'scrf',filesep,'movim2',ea_stripex(options.prefs.prenii_unnormalized),'_ants1.mat'],[directory,'scrf',filesep,'scrf_instore.mat']);
delete([directory,'scrf',filesep,ea_stripex(options.prefs.prenii_unnormalized),'2movim','_ants1.mat']);
ea_refreshscrf(options,handles,directory);




function otherfiles=ea_createmovim(directory,options)

switch options.modality
    case 1
        otherfiles={[directory,'scrf',filesep,options.prefs.tranii_unnormalized],...
            [directory,'scrf',filesep,options.prefs.cornii_unnormalized],...
            [directory,'scrf',filesep,options.prefs.sagnii_unnormalized]};
        cnt=1;
        for ofi=1:length(otherfiles)
            if exist(otherfiles{ofi},'file')
                nii=ea_load_nii(otherfiles{ofi});
                delete(otherfiles{ofi});
                nii.img(abs(nii.img)<0.1)=nan;
                if ~exist('AllX','var')
                    AllX=nii.img;
                else
                    AllX(:,:,:,cnt)=nii.img;
                end
                cnt=cnt+1;
            end
        end
        nii.img=ea_nanmean(AllX,4);
        clear AllX
        nii.fname=[directory,'scrf',filesep,'movim.nii'];
        nii.img(isnan(nii.img))=0;
        %nii.img(~(nii.img==0))=zscore(nii.img(~(nii.img==0)));

        ea_write_nii(nii);
    case 2
        otherfiles={[directory,'scrf',filesep,'tp_',options.prefs.ctnii_coregistered]};
        copyfile([directory,'scrf',filesep,'tp_',options.prefs.ctnii_coregistered],[directory,'scrf',filesep,'movim.nii']);
end

% --- Executes on button press in savebutn.
function savebutn_Callback(hObject, eventdata, handles)
% hObject    handle to savebutn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
directory=getappdata(handles.scrf,'directory');
if ~exist([directory,'scrf',filesep,'scrf_instore.mat'],'file')
	msgbox('Please generate a transform first (Click on "Compute subcortical refine transform"). If you don''t want to compute a transform, simply click on "Continue without subcortical transform".');
else
copyfile([directory,'scrf',filesep,'scrf_instore.mat'],[directory,'scrf',filesep,'scrf.mat']);
if exist([directory,'ea_reconstruction.mat'],'file')
ea_recalc_reco([],[],directory);
end

ea_methods(directory,...
            ['DBS electrode localizations were corrected for brainshift in postoperative acquisitions by applying a refined affine transform calculated between ',...
            'pre- and postoperative acquisitions that were restricted to a subcortical area of interest as implemented in the brainshift-correction module of Lead-DBS software',...
            ' (Horn & Kuehn 2005; SCR_002915; http://www.lead-dbs.org).'],...
            {'Horn, A., & Kühn, A. A. (2015). Lead-DBS: a toolbox for deep brain stimulation electrode localizations and visualizations. NeuroImage, 107, 127?135. http://doi.org/10.1016/j.neuroimage.2014.12.002'});

closescrf(handles);
end

% --- Executes on button press in rejectbutn.
function rejectbutn_Callback(hObject, eventdata, handles)
% hObject    handle to rejectbutn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
directory=getappdata(handles.scrf,'directory');
if exist([directory,'scrf',filesep,'scrf.mat']);
    delete([directory,'scrf',filesep,'scrf.mat']);
end
closescrf(handles);

function closescrf(handles)

options=getappdata(handles.scrf,'options');
if exist([options.root,options.patientname,filesep,'ea_reconstruction.mat'],'file') % apply brainshift correction to reconstruction
% read / write reco to include subcortical refine transform.

options.native=1;
[coords_mm,trajectory,markers,elmodel,manually_corrected]=ea_load_reconstruction(options);
options.hybridsave=1;
ea_save_reconstruction(coords_mm,trajectory,markers,elmodel,manually_corrected,options);
end
if options.d2.write || options.d3.write
uiresume(handles.scrf);
end
delete(handles.scrf);


% --- Executes on button press in mask1.
function mask1_Callback(hObject, eventdata, handles)
% hObject    handle to mask1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.mask0.Value=0;
handles.mask2.Value=0;
% Hint: get(hObject,'Value') returns toggle state of mask1


% --- Executes on button press in mask2.
function mask2_Callback(hObject, eventdata, handles)
% hObject    handle to mask2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.mask0.Value=0;
handles.mask1.Value=0;
% Hint: get(hObject,'Value') returns toggle state of mask2


% --- Executes on button press in mask0.
function mask0_Callback(hObject, eventdata, handles)
% hObject    handle to mask0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of mask0

handles.mask1.Value=0;
handles.mask2.Value=0;
% --- Executes on button press in mask1.
