function varargout = ea_imageclassifier(varargin)
% EA_IMAGECLASSIFIER MATLAB code for ea_imageclassifier.fig
%      EA_IMAGECLASSIFIER, by itself, creates a new EA_IMAGECLASSIFIER or raises the existing
%      singleton*.
%
%      H = EA_IMAGECLASSIFIER returns the handle to a new EA_IMAGECLASSIFIER or the handle to
%      the existing singleton*.
%
%      EA_IMAGECLASSIFIER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in EA_IMAGECLASSIFIER.M with the given input arguments.
%
%      EA_IMAGECLASSIFIER('Property','Value',...) creates a new EA_IMAGECLASSIFIER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ea_imageclassifier_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ea_imageclassifier_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ea_imageclassifier

% Last Modified by GUIDE v2.5 04-Jul-2016 12:39:31

% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ea_imageclassifier_OpeningFcn, ...
                   'gui_OutputFcn',  @ea_imageclassifier_OutputFcn, ...
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


% --- Executes just before ea_imageclassifier is made visible.
function ea_imageclassifier_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ea_imageclassifier (see VARARGIN)

% Choose default command line output for ea_imageclassifier
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);


set(hObject,'name','Please specify image type');

tmpoutdir=fileparts(varargin{1});

setappdata(hObject,'dcfilename',varargin{1});
setappdata(hObject,'tmpoutdir',tmpoutdir);

[pth,dcfn]=fileparts(getappdata(handles.imclassf,'dcfilename'));
% check for already classified images:
for append=[0,2:15] % check for duplicates, too.
    [base_lead_fis,all_lead_fis]=ea_getbasefilenames('',append);
    if ismember([dcfn,'.nii'],base_lead_fis)
        finishandclose(handles,'');
        return
    end
    if ismember([dcfn,'.nii'],all_lead_fis)
        finishandclose(handles,'');
        return
    end
    
    
    if append
        append=num2str(append);
    else
        append='';
    end
    
    if strcmp([dcfn,'.nii'],['unknown',append,'.nii'])
        finishandclose(handles,'');
        return
    end
    if strcmp([dcfn,'.nii'],['corrupt_data',append,'.nii'])
        finishandclose(handles,'');
        return
    end
    if strcmp([dcfn,'.nii'],['auto_ignore',append,'.nii'])
        finishandclose(handles,'');
        return
    end
end

try
nii=ea_load_nii(getappdata(handles.imclassf,'dcfilename'));
catch
    finishandclose(handles,'corrupt_data');
    return
end



if any(nii.dim<20)
finishandclose(handles,'auto_ignore');
return
end

hdrtext=genhdrtext(nii);


nii.img=double(nii.img)/double(max(nii.img(:)));
try
set(0,'CurrentFigure',handles.imclassf);
catch
    keyboard
end

try
h=slice(double(squeeze(nii.img(:,:,:,1))),round(size(nii.img,1)/2),...
    round(size(nii.img,2)/2),...
    round(size(nii.img,3)/2));
set(h,'FaceColor','interp',...
	'EdgeColor','none',...
	'DiffuseStrength',.8)
catch
    
    keyboard
end
colormap gray
set(hObject, 'menubar', 'figure' )
set(hObject, 'toolbar', 'figure' )
set(hObject,'KeyPressFcn',{@ea_keystr,handles});

view(45,40)
axis equal
axis tight
axis off

set(handles.imghdrinfo,'String',hdrtext);


% UIWAIT makes ea_imageclassifier wait for user response (see UIRESUME)
%uiwait(handles.imclassf);




% --- Outputs from this function are returned to the command line.
function varargout = ea_imageclassifier_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
%uiwait(handles.imclassf);



% --- Executes on button press in trapush.
function trapush_Callback(hObject, eventdata, handles)
% hObject    handle to trapush (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
tmpoutdir=getappdata(handles.imclassf,'tmpoutdir');
[~,pt]=fileparts(tmpoutdir);
prefs=ea_prefs(pt);
finishandclose(handles,prefs.tranii_unnormalized);

% --- Executes on button press in pretrapush.
function pretrapush_Callback(hObject, eventdata, handles)
% hObject    handle to pretrapush (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
tmpoutdir=getappdata(handles.imclassf,'tmpoutdir');
[~,pt]=fileparts(tmpoutdir);
prefs=ea_prefs(pt);
finishandclose(handles,prefs.prenii_unnormalized);

% --- Executes on button press in corpush.
function corpush_Callback(hObject, eventdata, handles)
% hObject    handle to corpush (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
tmpoutdir=getappdata(handles.imclassf,'tmpoutdir');
[~,pt]=fileparts(tmpoutdir);
prefs=ea_prefs(pt);
finishandclose(handles,prefs.cornii_unnormalized);

% --- Executes on button press in sagpush.
function sagpush_Callback(hObject, eventdata, handles)
% hObject    handle to sagpush (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
tmpoutdir=getappdata(handles.imclassf,'tmpoutdir');
[~,pt]=fileparts(tmpoutdir);
prefs=ea_prefs(pt);
finishandclose(handles,prefs.sagnii_unnormalized);

% --- Executes on button press in ctpush.
function ctpush_Callback(hObject, eventdata, handles)
% hObject    handle to ctpush (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
tmpoutdir=getappdata(handles.imclassf,'tmpoutdir');

[~,pt]=fileparts(tmpoutdir);
prefs=ea_prefs(pt);
finishandclose(handles,prefs.rawctnii_unnormalized);


% --- Executes on button press in unknownpush.
function unknownpush_Callback(hObject, eventdata, handles)
% hObject    handle to unknownpush (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
finishandclose(handles,'unknown');


function ea_keystr(icfig,event,handles)
%    pause
%commnd=get (handles.imclassf, 'CurrentKey');

        %global current_imclass
%% get vars
eltog=getappdata(handles.imclassf,'eltog');
elplot=getappdata(handles.imclassf,'elplot');
coords_mm=getappdata(handles.imclassf,'coords_mm');




commnd=lower(event.Character);
switch commnd
    case 't'
        tmpoutdir=getappdata(handles.imclassf,'tmpoutdir');
[~,pt]=fileparts(tmpoutdir);
prefs=ea_prefs(pt);
finishandclose(handles,prefs.tranii_unnormalized);
    case 'p'
tmpoutdir=getappdata(handles.imclassf,'tmpoutdir');
[~,pt]=fileparts(tmpoutdir);
prefs=ea_prefs(pt);
finishandclose(handles,prefs.prenii_unnormalized);
        
    case 's'
        tmpoutdir=getappdata(handles.imclassf,'tmpoutdir');
[~,pt]=fileparts(tmpoutdir);
prefs=ea_prefs(pt);
finishandclose(handles,prefs.sagnii_unnormalized);
        
        
    case 'c'
        tmpoutdir=getappdata(handles.imclassf,'tmpoutdir');
[~,pt]=fileparts(tmpoutdir);
prefs=ea_prefs(pt);
finishandclose(handles,prefs.cornii_unnormalized);
        
    case 'a'
        tmpoutdir=getappdata(handles.imclassf,'tmpoutdir');
[~,pt]=fileparts(tmpoutdir);
prefs=ea_prefs(pt);
finishandclose(handles,prefs.rawctnii_unnormalized);
        
    case 'u'
        finishandclose(handles,'unknown');
        
    case 'd'
        tmpoutdir=getappdata(handles.imclassf,'tmpoutdir');
        
        [~,pt]=fileparts(tmpoutdir);
        prefs=ea_prefs(pt);
        finishandclose(handles,prefs.dti);
    case 'f'
        tmpoutdir=getappdata(handles.imclassf,'tmpoutdir');
        
        [~,pt]=fileparts(tmpoutdir);
        prefs=ea_prefs(pt);
        finishandclose(handles,prefs.rest_default);
    otherwise
        return
end



function finishandclose(handles,current_imclass)
if ~isempty(current_imclass)
    [~,current_imclass]=fileparts(current_imclass); % remove potential file extension
    tmpoutdir=getappdata(handles.imclassf,'tmpoutdir');
    
    [~,pt]=fileparts(tmpoutdir);
    prefs=ea_prefs(pt);
    
    append='';
    while exist([tmpoutdir,filesep,current_imclass,append,'.nii'],'file')
        if isempty(append)
            append='2';
        else
            app=str2double(append);
            append=num2str(app+1);
        end
    end
    movefile([getappdata(handles.imclassf,'dcfilename')],[getappdata(handles.imclassf,'tmpoutdir'),filesep,current_imclass,append,'.nii']);
    
    [~,dti]=fileparts(prefs.dti);
    if strcmp(dti,current_imclass)
        
        fufn=getappdata(handles.imclassf,'dcfilename');
        [pth,fn,ext]=fileparts(fufn);
        if exist([pth,filesep,fn,'.bval'],'file')
            movefile([pth,filesep,fn,'.bval'],[pth,filesep,prefs.bval]);
        else
            warning('No .bval file found for dMRI image.');
        end
        
        if exist([pth,filesep,fn,'.bvec'],'file')
            movefile([pth,filesep,fn,'.bvec'],[pth,filesep,prefs.bvec]);
        else
            warning('No .bvec file found for dMRI image.');
        end
        
    end
    
end
close(handles.imclassf)
        


% --- Executes on button press in restpush.
function restpush_Callback(hObject, eventdata, handles)
% hObject    handle to restpush (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
tmpoutdir=getappdata(handles.imclassf,'tmpoutdir');

[~,pt]=fileparts(tmpoutdir);
prefs=ea_prefs(pt);
finishandclose(handles,prefs.rest_default);


% --- Executes on button press in dtipush.
function dtipush_Callback(hObject, eventdata, handles)
% hObject    handle to dtipush (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
tmpoutdir=getappdata(handles.imclassf,'tmpoutdir');

[~,pt]=fileparts(tmpoutdir);
prefs=ea_prefs(pt);
finishandclose(handles,prefs.dti);



function hdrtxt=genhdrtext(nii)

[~,fn,ext]=fileparts(nii.fname);

hdrtxt = '';
hdrtxt = sprintf('%sFile: \n', hdrtxt);
hdrtxt = sprintf('%s%s\n', hdrtxt, [fn,ext]);
hdrtxt = sprintf('%s----------------------------------\n', hdrtxt);
hdrtxt = sprintf('%sDimension: %d x %d x %d\n', hdrtxt, nii.dim(1),nii.dim(2),nii.dim(3));
hdrtxt = sprintf('%sIntensity Range: %d - %d\n', hdrtxt, ea_nanmin(nii.img(:)),ea_nanmax(nii.img(:)));
hdrtxt = sprintf('%sNumber of Components: %d\n', hdrtxt, size(nii.img,4));
hdrtxt = sprintf('%sVox2mm:\n', hdrtxt);
hdrtxt = sprintf('%s %.2f %.2f %.2f %.2f\n', hdrtxt,nii.mat(1,1),nii.mat(1,2),nii.mat(1,3),nii.mat(1,4));
hdrtxt = sprintf('%s %.2f %.2f %.2f %.2f\n', hdrtxt,nii.mat(2,1),nii.mat(2,2),nii.mat(2,3),nii.mat(2,4));
hdrtxt = sprintf('%s %.2f %.2f %.2f %.2f\n', hdrtxt,nii.mat(3,1),nii.mat(3,2),nii.mat(3,3),nii.mat(3,4));
hdrtxt = sprintf('%s %.2f %.2f %.2f %.2f\n', hdrtxt,nii.mat(4,1),nii.mat(4,2),nii.mat(4,3),nii.mat(4,4));


% --- Executes on button press in pretrat1push.
function pretrat1push_Callback(hObject, eventdata, handles)
% hObject    handle to pretrat1push (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
tmpoutdir=getappdata(handles.imclassf,'tmpoutdir');

[~,pt]=fileparts(tmpoutdir);
prefs=ea_prefs(pt);
finishandclose(handles,prefs.prenii_unnormalized_t1);

% --- Executes on button press in pretrapdpush.
function pretrapdpush_Callback(hObject, eventdata, handles)
% hObject    handle to pretrapdpush (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

tmpoutdir=getappdata(handles.imclassf,'tmpoutdir');

[~,pt]=fileparts(tmpoutdir);
prefs=ea_prefs(pt);
finishandclose(handles,prefs.prenii_unnormalized_pd);
