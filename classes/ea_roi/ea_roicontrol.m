function varargout = ea_roicontrol(varargin)
% EA_ROICONTROL MATLAB code for ea_roicontrol.fig
%      EA_ROICONTROL, by itself, creates a new EA_ROICONTROL or raises the existing
%      singleton*.
%
%      H = EA_ROICONTROL returns the handle to a new EA_ROICONTROL or the handle to
%      the existing singleton*.
%
%      EA_ROICONTROL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in EA_ROICONTROL.M with the given input arguments.
%
%      EA_ROICONTROL('Property','Value',...) creates a new EA_ROICONTROL or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ea_roicontrol_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ea_roicontrol_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ea_roicontrol

% Last Modified by GUIDE v2.5 25-Jun-2021 13:53:02

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @ea_roicontrol_OpeningFcn, ...
    'gui_OutputFcn',  @ea_roicontrol_OutputFcn, ...
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


% --- Executes just before ea_roicontrol is made visible.
function ea_roicontrol_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ea_roicontrol (see VARARGIN)

% Choose default command line output for ea_roicontrol
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ea_roicontrol wait for user response (see UIRESUME)
% uiwait(handles.roicontrol);
setappdata(handles.roicontrol,'chandles',handles);
obj=varargin{1};
setappdata(handles.roicontrol,'obj',obj);
nzeros=obj.nii.img(:);
nzeros(nzeros==0)=[];
nzeros(isnan(nzeros))=[];
set(0,'CurrentFigure',handles.roicontrol);
set(handles.roicontrol,'CurrentAxes',handles.histax);
axis off

handles.showhide.Value=ea_bool2onoff(obj.Visible);
handles.solidcolor.Visible=ea_bool2onoff(~obj.binary);
    
if ~obj.binary
    hist(nzeros);
    h=findobj(handles.histax,'Type','patch');
    set(h,'FaceColor',[0,0.5 0.5],'EdgeColor','none');
    axis off

    set(handles.roicontrol,'Name',obj.niftiFilename);
    handles.histax.XAxis.Limits=[min(nzeros),max(nzeros)];
else
    set(handles.threshLabel,'Visible','off');
end
% button
if ischar(obj.color) % none
    set(handles.colorchange,'BackgroundColor',[1,1,1]);
    setappdata(handles.colorchange,'facecolor',[1,1,1]);
    handles.facecolor.Value=0;
else
    set(handles.colorchange,'BackgroundColor',obj.color);
    setappdata(handles.colorchange,'facecolor',obj.color);
    handles.facecolor.Value=1;
end
if ischar(obj.edgecolor) % none
    set(handles.edgecolorchange,'BackgroundColor',[1,1,1]);
    setappdata(handles.edgecolorchange,'edgecolor',[1,1,1]);
    handles.edgecolor.Value=0;
else
    set(handles.edgecolorchange,'BackgroundColor',obj.edgecolor);
    setappdata(handles.edgecolorchange,'edgecolor',obj.edgecolor);
    handles.edgecolor.Value=1;
end
set(0,'CurrentFigure',handles.roicontrol);
set(handles.roicontrol,'name',obj.name);
%% sliders:
if ~obj.binary
    % threshold
    jSlider{1} = javax.swing.JSlider(0,100);
    ea_javacomponent(jSlider{1},[0,130,handles.roicontrol.Position(3)-1,45]);
    set(jSlider{1}, 'Value', getmaxminthresh(obj), 'MajorTickSpacing',0.1, 'PaintLabels',true);  % with labels, no ticks
    hjSlider{1} = handle(jSlider{1}, 'CallbackProperties');
    hjSlider{1}.Background=javax.swing.plaf.ColorUIResource(1,1,1);
    set(hjSlider{1}, 'MouseReleasedCallback', {@sliderthresholdchange,obj,handles});  %alternative
    set(hjSlider{1}, 'StateChangedCallback', {@sliderthresholdchangetxt,obj,handles});  %alternative
else
    obj.threshold=obj.max/2;
end
set(handles.threshtxt, 'String', num2str(obj.threshold));

% alpha
set(0,'CurrentFigure',handles.roicontrol);
jSlider{2} = javax.swing.JSlider(0,100);
ea_javacomponent(jSlider{2},[0,65,handles.roicontrol.Position(3)-1,45]);
set(jSlider{2}, 'Value', obj.alpha*100, 'MajorTickSpacing',0.1, 'PaintLabels',true);  % with labels, no ticks
set(handles.alphatxt, 'String', num2str(obj.alpha));
hjSlider{2} = handle(jSlider{2}, 'CallbackProperties');
hjSlider{2}.Background=javax.swing.plaf.ColorUIResource(1,1,1);
set(hjSlider{2}, 'StateChangedCallback', {@slideralphachange,obj,handles});  %alternative

% smooth
set(0,'CurrentFigure',handles.roicontrol);
jSlider{3} = javax.swing.JSlider(0,100);
ea_javacomponent(jSlider{3},[0,0,handles.roicontrol.Position(3)-1,45]);
set(jSlider{3}, 'Value', round(obj.smooth*2), 'MajorTickSpacing',0.1, 'PaintLabels',true);  % with labels, no ticks
set(handles.smoothtxt, 'String', [num2str(obj.smooth),' Its.']);
hjSlider{3} = handle(jSlider{3}, 'CallbackProperties');
hjSlider{3}.Background=javax.swing.plaf.ColorUIResource(1,1,1);
set(hjSlider{3}, 'StateChangedCallback', {@slidersmoothchangetxt,obj,handles});  %alternative
set(hjSlider{3}, 'MouseReleasedCallback', {@slidersmoothchange,obj,handles});  %alternative
set(0,'CurrentFigure',obj.plotFigureH);


function thresh=getmaxminthresh(obj)

thresh=obj.threshold-obj.min;
thresh=thresh/obj.max;
thresh=round(thresh*100); % slider only supports integers

function slideralphachange(varargin)
slide=varargin{1};
obj=varargin{3};
handles=varargin{4};
obj.alpha=slide.Value/100;
set(handles.alphatxt,'String',num2str(obj.alpha));

function sliderthresholdchange(varargin)
slide=varargin{1};
obj=varargin{3};
tval=slide.Value;
tval=tval/100;
tval=tval*obj.max;
obj.threshold=tval+obj.min;

function sliderthresholdchangetxt(varargin)
slide=varargin{1};
obj=varargin{3};
   handles=varargin{4};
tval=slide.Value;
tval=tval/100;
tval=tval*obj.max;
tval=tval+obj.min;
set(handles.threshtxt,'String',num2str(tval));

function slidersmoothchange(varargin)
slide=varargin{1};

obj=varargin{3};
obj.smooth=(round(slide.Value/2));

function slidersmoothchangetxt(varargin)
slide=varargin{1};
handles=varargin{4};
smooth=(round(slide.Value/2));
set(handles.smoothtxt,'String',[num2str(smooth),' Its.']);


% --- Outputs from this function are returned to the command line.
function varargout = ea_roicontrol_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in colorchange.
function colorchange_Callback(hObject, eventdata, handles)
% hObject    handle to colorchange (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
obj=getappdata(handles.roicontrol,'obj');
col = ea_uisetcolor;
setappdata(handles.colorchange,'facecolor',col);
set(handles.colorchange,'BackgroundColor',col);
if handles.facecolor.Value
    obj.color=getappdata(handles.colorchange,'facecolor');
end


% --- Executes on button press in showhide.
function showhide_Callback(hObject, eventdata, handles)
% hObject    handle to showhide (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of showhide
obj=getappdata(handles.roicontrol,'obj');
switch get(hObject,'Value')
    case 1
obj.Visible='on';
    case 0
   obj.Visible='off';
end


% --- Executes on button press in solidcolor.
function solidcolor_Callback(hObject, eventdata, handles)
% hObject    handle to solidcolor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of solidcolor
obj=getappdata(handles.roicontrol,'obj');
obj.usesolidcolor=get(hObject,'Value');


% --- Executes on button press in facecolor.
function facecolor_Callback(hObject, eventdata, handles)
% hObject    handle to facecolor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of facecolor
obj=getappdata(handles.roicontrol,'obj');
if handles.facecolor.Value
    obj.color=getappdata(handles.colorchange,'facecolor');
else
   obj.color='none';
end

% --- Executes on button press in edgecolorchange.
function edgecolorchange_Callback(hObject, eventdata, handles)
% hObject    handle to edgecolorchange (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
obj=getappdata(handles.roicontrol,'obj');
col = ea_uisetcolor;
setappdata(handles.edgecolorchange,'edgecolor',col);
set(handles.edgecolorchange,'BackgroundColor',col);
if handles.edgecolor.Value
    obj.edgecolor=getappdata(handles.edgecolorchange,'edgecolor');
end


% --- Executes on button press in edgecolor.
function edgecolor_Callback(hObject, eventdata, handles)
% hObject    handle to edgecolor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of edgecolor
obj=getappdata(handles.roicontrol,'obj');
if handles.edgecolor.Value
    obj.edgecolor=getappdata(handles.edgecolorchange,'edgecolor');
else
    obj.edgecolor='none';
end
