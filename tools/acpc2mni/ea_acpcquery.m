function varargout = ea_acpcquery(varargin)
% EA_ACPCQUERY MATLAB code for ea_acpcquery.fig
%      EA_ACPCQUERY, by itself, creates a new EA_ACPCQUERY or raises the existing
%      singleton*.
%
%      H = EA_ACPCQUERY returns the handle to a new EA_ACPCQUERY or the handle to
%      the existing singleton*.
%
%      EA_ACPCQUERY('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in EA_ACPCQUERY.M with the given input arguments.
%
%      EA_ACPCQUERY('Property','Value',...) creates a new EA_ACPCQUERY or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ea_acpcquery_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ea_acpcquery_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ea_acpcquery

% Last Modified by GUIDE v2.5 28-Jan-2016 11:55:48

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ea_acpcquery_OpeningFcn, ...
                   'gui_OutputFcn',  @ea_acpcquery_OutputFcn, ...
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


% --- Executes just before ea_acpcquery is made visible.
function ea_acpcquery_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ea_acpcquery (see VARARGIN)


earoot=ea_getearoot;
im=imread([earoot,'icons',filesep,'logo_lead_dbs.png']);
image(im);
axis off;
axis equal;

% Choose default command line output for ea_acpcquery
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ea_acpcquery wait for user response (see UIRESUME)
% uiwait(handles.acpcfig);
% UIWAIT makes tmp wait for user response (see UIRESUME)

setappdata(hObject,'leadfigure',varargin{3});
set(hObject,'name','ACPC/MNI-space conversions');

% --- Outputs from this function are returned to the command line.
function varargout = ea_acpcquery_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure



function xmm_Callback(hObject, eventdata, handles)
% hObject    handle to xmm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of xmm as text
%        str2double(get(hObject,'String')) returns contents of xmm as a double


% --- Executes during object creation, after setting all properties.
function xmm_CreateFcn(hObject, eventdata, handles)
% hObject    handle to xmm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ymm_Callback(hObject, eventdata, handles)
% hObject    handle to ymm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ymm as text
%        str2double(get(hObject,'String')) returns contents of ymm as a double


% --- Executes during object creation, after setting all properties.
function ymm_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ymm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function zmm_Callback(hObject, eventdata, handles)
% hObject    handle to zmm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of zmm as text
%        str2double(get(hObject,'String')) returns contents of zmm as a double


% --- Executes during object creation, after setting all properties.
function zmm_CreateFcn(hObject, eventdata, handles)
% hObject    handle to zmm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in xflip.
function xflip_Callback(hObject, eventdata, handles)
% hObject    handle to xflip (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns xflip contents as cell array
%        contents{get(hObject,'Value')} returns selected item from xflip


% --- Executes during object creation, after setting all properties.
function xflip_CreateFcn(hObject, eventdata, handles)
% hObject    handle to xflip (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in yflip.
function yflip_Callback(hObject, eventdata, handles)
% hObject    handle to yflip (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns yflip contents as cell array
%        contents{get(hObject,'Value')} returns selected item from yflip


% --- Executes during object creation, after setting all properties.
function yflip_CreateFcn(hObject, eventdata, handles)
% hObject    handle to yflip (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in zflip.
function zflip_Callback(hObject, eventdata, handles)
% hObject    handle to zflip (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns zflip contents as cell array
%        contents{get(hObject,'Value')} returns selected item from zflip


% --- Executes during object creation, after setting all properties.
function zflip_CreateFcn(hObject, eventdata, handles)
% hObject    handle to zflip (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in ac.
function ac_Callback(hObject, eventdata, handles)
% hObject    handle to ac (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ac
set(handles.mc,'Value',~get(handles.ac,'Value'));
set(handles.pc,'Value',~get(handles.ac,'Value'));

% --- Executes on button press in mc.
function mc_Callback(hObject, eventdata, handles)
% hObject    handle to mc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of mc
set(handles.ac,'Value',~get(handles.mc,'Value'));
set(handles.pc,'Value',~get(handles.mc,'Value'));

% --- Executes on button press in pc.
function pc_Callback(hObject, eventdata, handles)
% hObject    handle to pc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.ac,'Value',~get(handles.pc,'Value'));
set(handles.mc,'Value',~get(handles.pc,'Value'));
% Hint: get(hObject,'Value') returns toggle state of pc


% --- Executes on button press in cancelbutn.
function cancelbutn_Callback(hObject, eventdata, handles)
% hObject    handle to cancelbutn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.output = 'canceled';

% Update handles structure
guidata(hObject, handles);

% Use UIRESUME instead of delete because the OutputFcn needs
% to get the updated handles structure.
uiresume(handles.acpcfig);

% --- Executes on button press in acpc2mnibutn.
function acpc2mnibutn_Callback(hObject, eventdata, handles)
% hObject    handle to acpc2mnibutn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

ea_busyaction('on',handles.acpcfig,'acpc');

cfg.xmm=str2double(get(handles.xmm,'String'));
cfg.ymm=str2double(get(handles.ymm,'String'));
cfg.zmm=str2double(get(handles.zmm,'String'));
if get(handles.xflip,'Value')==2
    cfg.xmm=cfg.xmm*-1;
end
if get(handles.yflip,'Value')==2
    cfg.ymm=cfg.ymm*-1;
end
if get(handles.zflip,'Value')==1
    cfg.zmm=cfg.zmm*-1;
end
if get(handles.ac,'Value')
    cfg.acmcpc=1;
elseif get(handles.mc,'Value')
    cfg.acmcpc=2;
elseif get(handles.pc,'Value')
    cfg.acmcpc=3;
end

cfg.mapmethod=get(handles.methodm,'Value')-1;

handles.output = cfg;


% Update handles structure
guidata(hObject, handles);

% Use UIRESUME instead of delete because the OutputFcn needs
% to get the updated handles structure.
leadfigure=getappdata(handles.acpcfig,'leadfigure');

fid=ea_acpc2mni(cfg,leadfigure);


leaddir=[ea_getearoot];
        tempfile=[leaddir,'templates',filesep,'mni_hires.nii'];


for pt=1:length(fid)
    mnipoints(pt,:)=fid(pt).WarpedPointMNI;
end

meanmni=mean(mnipoints,1);
set(handles.xmni,'String',num2str(meanmni(1))); set(handles.ymni,'String',num2str(meanmni(2))); set(handles.zmni,'String',num2str(meanmni(3)));
stdmni=std(mnipoints,0,1);
set(handles.xstdmni,'String',['± ',sprintf('%.2f', stdmni(1)),' mm']); set(handles.ystdmni,'String',['± ',sprintf('%.2f', stdmni(2)),' mm']); set(handles.zstdmni,'String',['± ',sprintf('%.2f', stdmni(3)),' mm']);

ea_busyaction('off',handles.acpcfig,'acpc');


% --- Executes on selection change in methodm.
function methodm_Callback(hObject, eventdata, handles)
% hObject    handle to methodm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns methodm contents as cell array
%        contents{get(hObject,'Value')} returns selected item from methodm


% --- Executes during object creation, after setting all properties.
function methodm_CreateFcn(hObject, eventdata, handles)
% hObject    handle to methodm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function xmni_Callback(hObject, eventdata, handles)
% hObject    handle to xmni (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of xmni as text
%        str2double(get(hObject,'String')) returns contents of xmni as a double


% --- Executes during object creation, after setting all properties.
function xmni_CreateFcn(hObject, eventdata, handles)
% hObject    handle to xmni (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ymni_Callback(hObject, eventdata, handles)
% hObject    handle to ymni (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ymni as text
%        str2double(get(hObject,'String')) returns contents of ymni as a double


% --- Executes during object creation, after setting all properties.
function ymni_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ymni (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function zmni_Callback(hObject, eventdata, handles)
% hObject    handle to zmni (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of zmni as text
%        str2double(get(hObject,'String')) returns contents of zmni as a double


% --- Executes during object creation, after setting all properties.
function zmni_CreateFcn(hObject, eventdata, handles)
% hObject    handle to zmni (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in mni2acpcbutn.
function mni2acpcbutn_Callback(hObject, eventdata, handles)
% hObject    handle to mni2acpcbutn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


ea_busyaction('on',handles.acpcfig,'acpc');

cfg.xmm=str2double(get(handles.xmni,'String'));
cfg.ymm=str2double(get(handles.ymni,'String'));
cfg.zmm=str2double(get(handles.zmni,'String'));

if get(handles.ac,'Value')
    cfg.acmcpc=1;
elseif get(handles.mc,'Value')
    cfg.acmcpc=2;
elseif get(handles.pc,'Value')
    cfg.acmcpc=3;
end



% Use UIRESUME instead of delete because the OutputFcn needs
% to get the updated handles structure.
leadfigure=getappdata(handles.acpcfig,'leadfigure');
fid=ea_mni2acpc(cfg,leadfigure);






for pt=1:length(fid)
    acpcpoints(pt,:)=fid(pt).WarpedPointACPC;
end

meanacpc=mean(acpcpoints,1);

stdacpc=std(acpcpoints,0,1);

% default: 112

if get(handles.xflip,'Value')==2
    meanacpc(1)=meanacpc(1)*-1;
end
if get(handles.yflip,'Value')==2
    meanacpc(2)=meanacpc(2)*-1;
end
if get(handles.zflip,'Value')==1
    meanacpc(3)=meanacpc(3)*-1;
end

set(handles.xmm,'String',num2str(meanacpc(1))); set(handles.ymm,'String',num2str(meanacpc(2))); set(handles.zmm,'String',num2str(meanacpc(3)));
set(handles.xstdacpc,'String',['± ',sprintf('%.2f', stdacpc(1)),' mm']); set(handles.ystdacpc,'String',['± ',sprintf('%.2f', stdacpc(2)),' mm']); set(handles.zstdacpc,'String',['± ',sprintf('%.2f', stdacpc(3)),' mm']);

ea_busyaction('off',handles.acpcfig,'acpc');
