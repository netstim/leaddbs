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

% Last Modified by GUIDE v2.5 03-Mar-2017 15:41:22

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
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




directory=[varargin{1},filesep];
lfhandles=varargin{2};
setappdata(handles.scrf,'lfhandles',lfhandles);
setappdata(handles.scrf,'directory',directory);

[pth,patientname]=fileparts(fileparts(directory));
handles.patientname.String=patientname;

ea_refreshscrf(lfhandles,handles,directory);


% Choose default command line output for ea_subcorticalrefine
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ea_subcorticalrefine wait for user response (see UIRESUME)
% uiwait(handles.scrf);

function ea_refreshscrf(lfhandles,handles,directory)

standardslice=ea_loadrefineslice(directory,lfhandles,0);
refineslice=ea_loadrefineslice(directory,lfhandles,1);
set(0,'CurrentFigure',handles.scrf);

handles.scrf.CurrentAxes=handles.standardax;
imshow(standardslice);
handles.scrf.CurrentAxes=handles.scfax;
imshow(refineslice);

function slice=ea_loadrefineslice(directory,lfhandles,refine)

switch refine
    case 1
        refstr='scrf';
    case 0
        refstr='standard';
end

ea_createrefineslice(directory,lfhandles,refine); 
%ea_createcompimg(directory,refine,refstr,lfhandles)
try
slice=imread([directory,'scrf',filesep,refstr,'.png']);
catch
    slice=zeros(10,10,3);
end

function ea_createcompimg(directory,refine,refstr,lfhandles)
return
options.modality=lfhandles.MRCT.Value;
[options.root,options.patientname]=fileparts(fileparts(directory));
options.root=[options.root,filesep];
options.prefs=ea_prefs(options.patientname);
options=ea_assignpretra(options);
switch refine
    case 1
        scrf='scrf';
    case 0
        scrf='';
end
switch options.modality
    case 1
        fis={options.prefs.tranii_unnormalized,options.prefs.cornii_unnormalized,options.prefs.sagnii_unnormalized};
        
        for fi=1:length(fis)
            cnt=1;
            if exist([directory,'scrf',filesep,scrf,ea_stripex(fis{fi}),'2',ea_stripex(options.prefs.prenii_unnormalized),'.png'],'file');
                thisim=imread([directory,'scrf',filesep,scrf,ea_stripex(fis{fi}),'2',ea_stripex(options.prefs.prenii_unnormalized),'.png']);
                if ~exist('chan','var')
                    chan=thisim;
                else
                    chan=chan+thisim;
                end
                cnt=cnt+1;
            end
            if ~exist('chan','var') % files not present.
                chan=zeros(100,100,3);
            end
            chan=chan/cnt;
        end
    case 2
        chan=imread([directory,'scrf',filesep,scrf,'tp_',ea_stripex(options.prefs.ct_coregistered),'2',ea_stripex(options.prefs.prenii_unnormalized),'.png']);

end

imwrite(chan,[directory,'scrf',filesep,refstr,'.png']);


function fn=ea_stripex(fn)
[~,fn]=fileparts(fn);

function ea_createrefineslice(directory,lfhandles,refine)


options.modality=lfhandles.MRCT.Value;
[options.root,options.patientname]=fileparts(fileparts(directory));
options.root=[options.root,filesep];
options.prefs=ea_prefs(options.patientname);
options=ea_assignpretra(options);
switch refine
    case 1
        scrf='scrf';
    case 0
        scrf='';
end


if ~exist([directory,'scrf',filesep,options.prefs.prenii_unnormalized],'file')
    ea_createbbfiles(directory);
end
ea_createmovim(directory,options);
ea_gencoregcheckfigs_scrf(directory,scrf,options);


function ea_createbbfiles(directory)

[options.root,options.patientname]=fileparts(directory);
options.root=[options.root,filesep];
options.earoot=ea_getearoot;
options.prefs=ea_prefs(options.patientname);
options=ea_assignpretra(options);

if ~exist([directory,'scrf'],'dir')
    mkdir([directory,'scrf']);
end

to{1}=[directory,'scrf',filesep,'bb.nii'];
from{1}=[ea_space,'bb.nii'];


ea_apply_normalization_tofile(options,from,to,[options.root,options.patientname,filesep],1);
ea_crop_nii([directory,'scrf',filesep,'bb.nii']);
fis={options.prefs.prenii_unnormalized,options.prefs.tranii_unnormalized,options.prefs.cornii_unnormalized,options.prefs.sagnii_unnormalized,options.prefs.ctnii_coregistered};
for fi=1:length(fis)
    if exist([directory,fis{fi}],'file')
        copyfile([directory,fis{fi}],[directory,'scrf',filesep,fis{fi}])
        ea_conformspaceto([directory,'scrf',filesep,'bb.nii'],[directory,'scrf',filesep,fis{fi}]);
    end
end

% --- Outputs from this function are returned to the command line.
function varargout = ea_subcorticalrefine_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


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

lfhandles=getappdata(handles.scrf,'lfhandles');
options.modality=lfhandles.MRCT.Value;
directory=getappdata(handles.scrf,'directory');

[options.root,options.patientname]=fileparts(fileparts(directory));
options.root=[options.root,filesep];
   options.prefs=ea_prefs('');
   options=ea_assignpretra(options);
   
if iscell(handles.methodbutton.String)
    method=handles.methodbutton.String{handles.methodbutton.Value};
else
    method=handles.methodbutton.String;
end
switch(method)
    case 'ANTs'
     options.coregmr.method='Coreg MRIs: ANTs';
    case 'SPM'
     options.coregmr.method='Coreg MRIs: SPM';  
end
otherfiles=ea_createmovim(directory,options);


ea_coreg2images(options,[directory,'scrf',filesep,'movim.nii'],[directory,'scrf',filesep,options.prefs.prenii_unnormalized],...
    [directory,'scrf',filesep,'scrfmovim.nii'],otherfiles,1);
movefile([directory,'scrf',filesep,'movim.nii'],[directory,'scrf',filesep,'scrfmovim.nii']);
movefile([directory,'scrf',filesep,'raw_movim.nii'],[directory,'scrf',filesep,'movim.nii']);
ea_refreshscrf(lfhandles,handles,directory);

movefile([directory,'scrf',filesep,'movim2',ea_stripex(options.prefs.prenii_unnormalized),'_ants1.mat'],[directory,'scrf',filesep,'scrf_instore.mat']);
delete([directory,'scrf',filesep,ea_stripex(options.prefs.prenii_unnormalized),'2movim','_ants1.mat']);

function otherfiles=ea_createmovim(directory,options)

switch options.modality
    case 1
        otherfiles={[directory,'scrf',filesep,options.prefs.tranii_unnormalized],...
            [directory,'scrf',filesep,options.prefs.cornii_unnormalized],...
            [directory,'scrf',filesep,options.prefs.sagnii_unnormalized]};
        
        for ofi=1:length(otherfiles)
           nii=ea_load_nii(otherfiles{ofi});
           if ~exist('AllX','var')
               AllX=nii.img;
           else
               AllX=AllX+nii.img;
           end
        end
        nii.img=AllX./ofi;
        clear AllX
        nii.fname=[directory,'scrf',filesep,'movim.nii'];
        ea_write_nii(nii);
    case 2
        otherfiles={[directory,'scrf',filesep,options.prefs.ctnii_coregistered]};
        copyfile([directory,'scrf',filesep,options.prefs.ctnii_coregistered],[directory,'scrf',filesep,'movim.nii']);
end

% --- Executes on button press in savebutn.
function savebutn_Callback(hObject, eventdata, handles)
% hObject    handle to savebutn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
directory=getappdata(handles.scrf,'directory');
copyfile([directory,'scrf',filesep,'scrf_instore.mat'],[directory,'scrf',filesep,'scrf.mat']);
handles.msgtxt.String='Transform saved.';


% --- Executes on button press in rejectbutn.
function rejectbutn_Callback(hObject, eventdata, handles)
% hObject    handle to rejectbutn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
