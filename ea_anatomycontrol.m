function varargout = ea_anatomycontrol(varargin)
% EA_ANATOMYCONTROL MATLAB code for ea_anatomycontrol.fig
%      EA_ANATOMYCONTROL, by itself, creates a new EA_ANATOMYCONTROL or raises the existing
%      singleton*.
%
%      H = EA_ANATOMYCONTROL returns the handle to a new EA_ANATOMYCONTROL or the handle to
%      the existing singleton*.
%
%      EA_ANATOMYCONTROL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in EA_ANATOMYCONTROL.M with the given input arguments.
%
%      EA_ANATOMYCONTROL('Property','Value',...) creates a new EA_ANATOMYCONTROL or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ea_anatomycontrol_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ea_anatomycontrol_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ea_anatomycontrol

% Last Modified by GUIDE v2.5 08-Feb-2016 18:56:44

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ea_anatomycontrol_OpeningFcn, ...
                   'gui_OutputFcn',  @ea_anatomycontrol_OutputFcn, ...
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


% --- Executes just before ea_anatomycontrol is made visible.
function ea_anatomycontrol_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ea_anatomycontrol (see VARARGIN)

set(hObject,'Name','Anatomy Slices');

% Choose default command line output for ea_anatomycontrol
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ea_anatomycontrol wait for user response (see UIRESUME)
% uiwait(handles.acontrolfig);
resultfig=varargin{1};
options=varargin{2};
setappdata(hObject,'resultfig',resultfig);
setappdata(hObject,'options',options);
togglestates=getappdata(resultfig,'togglestates'); % get info from resultfig.
setappdata(hObject,'togglestates',togglestates); % store anatomy toggle data from resultfig to anatomyslice (this) fig for subroutines.
set(handles.acontrolfig,'Visible',options.d3.verbose);


if ~isfield(options,'native')
    options.native=0;
end
list=ea_assignbackdrop('list',options,'Patient',options.native);
set(handles.templatepopup,'String',list);

if ~isempty(togglestates) % anatomy toggles have been used before..
% reset figure handle.

% xyz values
set(handles.xval,'String',num2str(togglestates.xyzmm(1)));
set(handles.yval,'String',num2str(togglestates.xyzmm(2)));
set(handles.zval,'String',num2str(togglestates.xyzmm(3)));
% toggle buttons
set(handles.xtoggle,'Value',togglestates.xyztoggles(1));
set(handles.ytoggle,'Value',togglestates.xyztoggles(2));
set(handles.ztoggle,'Value',togglestates.xyztoggles(3));
% transparencies
set(handles.xtrans,'String',num2str(togglestates.xyztransparencies(1)));
set(handles.ytrans,'String',num2str(togglestates.xyztransparencies(2)));
set(handles.ztrans,'String',num2str(togglestates.xyztransparencies(3)));

% template name
set(handles.templatepopup,'Value',find(ismember(get(handles.templatepopup,'String'),togglestates.template)));

% invertcheck
set(handles.invertcheck,'Value',togglestates.tinvert);

% cut items.
switch togglestates.cutview
    case '3d'
        set(handles.threedradio,'Value',1);
    case 'xcut'
        set(handles.xcutradio,'Value',1);
    case 'ycut'
        set(handles.ycutradio,'Value',1);
    case 'zcut'
        set(handles.zcutradio,'Value',1);
end
else
    togglestates.cutview='3d';
    
setappdata(getappdata(handles.acontrolfig,'resultfig'),'togglestates',togglestates);
end

pos=get(hObject,'position');
set(hObject,'position',[0,0,pos(3),pos(4)]);
refreshresultfig(handles)






% --- Outputs from this function are returned to the command line.
function varargout = ea_anatomycontrol_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



% --- Executes on selection change in templatepopup.
function templatepopup_Callback(hObject, eventdata, handles)
% hObject    handle to templatepopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns templatepopup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from templatepopup
popvals=get(hObject,'String');
if strcmp(popvals{get(hObject,'Value')},'Choose...')
    [FileName,PathName,FilterIndex] = uigetfile('*.nii','Choose anatomical image...');
    setappdata(gcf,'customfile',[PathName,FileName]);
end
    
refreshresultfig(handles)

% --- Executes during object creation, after setting all properties.
function templatepopup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to templatepopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in xtoggle.
function xtoggle_Callback(hObject, eventdata, handles)
% hObject    handle to xtoggle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
refreshresultfig(handles)

% --- Executes on button press in ytoggle.
function ytoggle_Callback(hObject, eventdata, handles)
% hObject    handle to ytoggle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
refreshresultfig(handles)


% --- Executes on button press in ztoggle.
function ztoggle_Callback(hObject, eventdata, handles)
% hObject    handle to ztoggle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
refreshresultfig(handles)


function xval_Callback(hObject, eventdata, handles)
% hObject    handle to xval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of xval as text
%        str2double(get(hObject,'String')) returns contents of xval as a double

refreshresultfig(handles)

% --- Executes during object creation, after setting all properties.
function xval_CreateFcn(hObject, eventdata, handles)
% hObject    handle to xval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function yval_Callback(hObject, eventdata, handles)
% hObject    handle to yval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of yval as text
%        str2double(get(hObject,'String')) returns contents of yval as a double
refreshresultfig(handles)

% --- Executes during object creation, after setting all properties.
function yval_CreateFcn(hObject, eventdata, handles)
% hObject    handle to yval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function zval_Callback(hObject, eventdata, handles)
% hObject    handle to zval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of zval as text
%        str2double(get(hObject,'String')) returns contents of zval as a double
refreshresultfig(handles)

% --- Executes during object creation, after setting all properties.
function zval_CreateFcn(hObject, eventdata, handles)
% hObject    handle to zval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function xtrans_Callback(hObject, eventdata, handles)
% hObject    handle to xtrans (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of xtrans as text
%        str2double(get(hObject,'String')) returns contents of xtrans as a double
refreshresultfig(handles)


% --- Executes during object creation, after setting all properties.
function xtrans_CreateFcn(hObject, eventdata, handles)
% hObject    handle to xtrans (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ytrans_Callback(hObject, eventdata, handles)
% hObject    handle to ytrans (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ytrans as text
%        str2double(get(hObject,'String')) returns contents of ytrans as a double
refreshresultfig(handles)


% --- Executes during object creation, after setting all properties.
function ytrans_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ytrans (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ztrans_Callback(hObject, eventdata, handles)
% hObject    handle to ztrans (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ztrans as text
%        str2double(get(hObject,'String')) returns contents of ztrans as a double
refreshresultfig(handles)


% --- Executes during object creation, after setting all properties.
function ztrans_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ztrans (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in invertcheck.
function invertcheck_Callback(hObject, eventdata, handles)
% hObject    handle to invertcheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of invertcheck
refreshresultfig(handles)


function refreshresultfig(handles)
% this part makes changes of the figure active:


resultfig=getappdata(handles.acontrolfig,'resultfig');
try
togglestates=getappdata(resultfig,'togglestates');
catch
    resultfig=gcf;
togglestates=getappdata(handles.acontrolfig,'togglestates');
end
if ~isfield(togglestates,'refreshcuts')
    togglestates.refreshcuts=0;
end

% reset states based on gui:
togglestates.xyzmm=[str2double(get(handles.xval,'String')),str2double(get(handles.yval,'String')),str2double(get(handles.zval,'String'))];
togglestates.xyztoggles=[get(handles.xtoggle,'Value'),get(handles.ytoggle,'Value'),get(handles.ztoggle,'Value')];
togglestates.xyztransparencies=[str2double(get(handles.xtrans,'String')),str2double(get(handles.ytrans,'String')),str2double(get(handles.ztrans,'String'))];
togglestates.template=get(handles.templatepopup,'String');
togglestates.template=togglestates.template{get(handles.templatepopup,'Value')};
togglestates.tinvert=get(handles.invertcheck,'Value');
togglestates.customfile=getappdata(gcf,'customfile');
setappdata(getappdata(handles.acontrolfig,'resultfig'),'togglestates',togglestates); % also store toggle data in resultfig.

ea_anatomyslices(getappdata(handles.acontrolfig,'resultfig'),...
    togglestates,...
    getappdata(handles.acontrolfig,'options'));

switch togglestates.cutview
    case 'xcut'
        set(0,'CurrentFigure',resultfig);
                [az,el]=view;
        if ~(az==90 && el==0)  || togglestates.refreshcuts          
            axis([str2double(get(handles.xval,'String'))-1 str2double(get(handles.xval,'String'))+1 -130 100 -70 100])
            view(90,0);
        end
    case 'ycut'
        set(0,'CurrentFigure',resultfig);
                [az,el]=view;
        if ~(az==0 && el==0) || togglestates.refreshcuts
            axis([-100 100 str2double(get(handles.yval,'String'))-1 str2double(get(handles.yval,'String'))+1 -70 100])
            view(0,0);
        end
    case 'zcut'
        set(0,'CurrentFigure',resultfig);
        [az,el]=view;
        if ~(az==0 && el==90) || togglestates.refreshcuts
            axis([-100 100 -130 100 str2double(get(handles.zval,'String'))-1,str2double(get(handles.zval,'String'))+1])
            view(0,90);
        end
    case '3d'
        set(0,'CurrentFigure',resultfig);
        if togglestates.refreshcuts; axis([-100 100 -130 100 -70 100]); end
        
end

togglestates.refreshcuts=0;
setappdata(resultfig,'togglestates',togglestates);


% --------------------------------------------------------------------
function slicebuttongroup_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to slicebuttongroup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)




% --- Executes during object creation, after setting all properties.
function slicebuttongroup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slicebuttongroup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object deletion, before destroying properties.
function slicebuttongroup_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to slicebuttongroup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes when slicebuttongroup is resized.
function slicebuttongroup_ResizeFcn(hObject, eventdata, handles)
% hObject    handle to slicebuttongroup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes when selected object is changed in slicebuttongroup.
function slicebuttongroup_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in slicebuttongroup 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)


% always and only store togglestates in the *resultfig*, not in the small
% control figure (where only the handle to resultfig is stored.
togglestates=getappdata(getappdata(gcf,'resultfig'),'togglestates');

switch eventdata.NewValue
    case handles.xcutradio
        togglestates.cutview='xcut';
        if eventdata.OldValue==handles.threedradio; togglestates.refreshcuts=1; end
            
    case handles.ycutradio
        togglestates.cutview='ycut';
        if eventdata.OldValue==handles.threedradio; togglestates.refreshcuts=1; end
        
    case handles.zcutradio
        togglestates.cutview='zcut';
        if eventdata.OldValue==handles.threedradio; togglestates.refreshcuts=1; end
        
    case handles.threedradio
        
        togglestates.cutview='3d';
        togglestates.refreshcuts=1;
        
end
setappdata(getappdata(gcf,'resultfig'),'togglestates',togglestates);

refreshresultfig(handles)


% --- Executes on button press in specify2dwrite.
function specify2dwrite_Callback(hObject, eventdata, handles)
% hObject    handle to specify2dwrite (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ea_spec2dwrite;
