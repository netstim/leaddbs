function varargout = ea_sweetspotcontrol(varargin)
% EA_SWEETSPOTCONTROL MATLAB code for ea_sweetspotcontrol.fig
%      EA_SWEETSPOTCONTROL, by itself, creates a new EA_SWEETSPOTCONTROL or raises the existing
%      singleton*.
%
%      H = EA_SWEETSPOTCONTROL returns the handle to a new EA_SWEETSPOTCONTROL or the handle to
%      the existing singleton*.
%
%      EA_SWEETSPOTCONTROL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in EA_SWEETSPOTCONTROL.M with the given input arguments.
%
%      EA_SWEETSPOTCONTROL('Property','Value',...) creates a new EA_SWEETSPOTCONTROL or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ea_sweetspotcontrol_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ea_sweetspotcontrol_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ea_sweetspotcontrol

% Last Modified by GUIDE v2.5 07-Feb-2018 14:15:38

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @ea_sweetspotcontrol_OpeningFcn, ...
    'gui_OutputFcn',  @ea_sweetspotcontrol_OutputFcn, ...
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


% --- Executes just before ea_sweetspotcontrol is made visible.
function ea_sweetspotcontrol_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ea_sweetspotcontrol (see VARARGIN)

% Choose default command line output for ea_sweetspotcontrol
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ea_sweetspotcontrol wait for user response (see UIRESUME)
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

if ~obj.binary
    hist(nzeros);
    h=findobj(handles.histax,'Type','patch');
    set(h,'FaceColor',[0,0.5 0.5],'EdgeColor','none');
    axis off

    set(handles.roicontrol,'Name',obj.niftiFilename);
    handles.histax.XAxis.Limits=[min(nzeros),max(nzeros)];
else
    set(handles.threshtxt,'Visible','off');
end
% button
set(handles.colorchange,'BackgroundColor',obj.color);
set(0,'CurrentFigure',handles.roicontrol);
set(handles.roicontrol,'name',obj.name);
%% sliders:
if ~obj.binary
% threshold
jSlider{1} = javax.swing.JSlider(0,100);
javacomponent(jSlider{1},[0,130,200,45]);
set(jSlider{1}, 'Value', getmaxminthresh(obj), 'MajorTickSpacing',0.1, 'PaintLabels',true);  % with labels, no ticks
hjSlider{1} = handle(jSlider{1}, 'CallbackProperties');
set(hjSlider{1}, 'MouseReleasedCallback', {@sliderthresholdchange,obj,handles});  %alternative
set(hjSlider{1}, 'StateChangedCallback', {@sliderthresholdchangetxt,obj,handles});  %alternative
else
   obj.threshold=obj.max/2;
end

% alpha
set(0,'CurrentFigure',handles.roicontrol);
jSlider{2} = javax.swing.JSlider(0,100);
javacomponent(jSlider{2},[0,65,200,45]);
set(jSlider{2}, 'Value', obj.alpha*100, 'MajorTickSpacing',0.1, 'PaintLabels',true);  % with labels, no ticks
hjSlider{2} = handle(jSlider{2}, 'CallbackProperties');
set(hjSlider{2}, 'StateChangedCallback', {@slideralphachange,obj,handles});  %alternative

% smooth
set(0,'CurrentFigure',handles.roicontrol);
jSlider{3} = javax.swing.JSlider(0,100);
javacomponent(jSlider{3},[0,0,200,45]);
set(jSlider{3}, 'Value', round(obj.smooth*2), 'MajorTickSpacing',0.1, 'PaintLabels',true);  % with labels, no ticks
hjSlider{3} = handle(jSlider{3}, 'CallbackProperties');
set(hjSlider{3}, 'StateChangedCallback', {@slidersmoothchangetxt,obj,handles});  %alternative
set(hjSlider{3}, 'MouseReleasedCallback', {@slidersmoothchange,obj,handles});  %alternative
set(0,'CurrentFigure',obj.plotFigureH);

cnt=4+size(obj.nii.X,1);
% factors
fcnt=1;
for f=size(obj.nii.X,1):-1:1
   set(0,'CurrentFigure',handles.roicontrol);

            jSlider{cnt} = javax.swing.JSlider(0,100);
javacomponent(jSlider{cnt},[200,0+(fcnt-1)*65,200,45]);
set(jSlider{cnt}, 'Value', 100/length(obj.nii.factors), 'MajorTickSpacing',0.1, 'PaintLabels',true);  % with labels, no ticks
FactorSlider{f} = handle(jSlider{cnt}, 'CallbackProperties');
set(FactorSlider{f}, 'MouseReleasedCallback', {@sliderfactorchange,obj,handles,f});  %alternative
set(FactorSlider{f}, 'StateChangedCallback', {@sliderfactorchangelive,obj,handles,f});  %alternative   
t = uicontrol(handles.roicontrol,'Style','text',...
                'String',obj.nii.factors{f},...
                'Position',[200,(0+(fcnt-1)*65)+50,200,12]);
cnt=cnt-1;
fcnt=fcnt+1;
end
setappdata(handles.roicontrol,'FactorSlider',FactorSlider);
    

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

function sliderfactorchangelive(varargin)

slide=varargin{1};
obj=varargin{3};
handles=varargin{4};
factor=varargin{5};
FactorSlider=getappdata(handles.roicontrol,'FactorSlider');

for slider=1:length(FactorSlider)
tval(slider)=FactorSlider{slider}.Value;
tval(slider)=tval(slider)/100;
end
% update fmix model:

ofactors=1:length(obj.nii.factors); ofactors(factor)=[];
tval(ofactors)=(tval(ofactors)/sum(tval(ofactors)))*(1-tval(factor));
tval(isnan(tval))=eps;
% scale remaining sliders accordingly:
for osliders=ofactors
    FactorSlider{osliders}.Value=tval(osliders)*100;
end

function sliderfactorchange(varargin)

slide=varargin{1};
obj=varargin{3};
handles=varargin{4};
factor=varargin{5};
FactorSlider=getappdata(handles.roicontrol,'FactorSlider');

tval=slide.Value;
tval=tval/100;
% update fmix model:
obj.fmix(factor)=tval;
ofactors=1:length(obj.nii.factors); ofactors(factor)=[];
obj.fmix(ofactors)=(obj.fmix(ofactors)/sum(obj.fmix(ofactors)))*(1-tval);
obj.fmix(isnan(obj.fmix))=eps;
% scale remaining sliders accordingly:
for osliders=ofactors
    FactorSlider{osliders}.Value=obj.fmix(osliders)*100;
end


function sliderthresholdchangetxt(varargin)
slide=varargin{1};
obj=varargin{3};
   handles=varargin{4};
tval=slide.Value;
tval=tval/100;
tval=tval*obj.max;
tval=tval+obj.min;
set(handles.theshtxt,'String',sprintf('%0.2f',tval));

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
function varargout = ea_sweetspotcontrol_OutputFcn(hObject, eventdata, handles)
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
obj.color = ea_uisetcolor;
set(handles.colorchange,'BackgroundColor',obj.color);


% --- Executes on button press in showhide.
function showhide_Callback(hObject, eventdata, handles)
% hObject    handle to showhide (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of showhide
obj=getappdata(handles.roicontrol,'obj');
switch get(hObject,'Value')
    case 1
obj.visible='on';
    case 0
   obj.visible='off';
end
