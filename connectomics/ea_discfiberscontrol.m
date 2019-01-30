function varargout = ea_discfiberscontrol(varargin)
% EA_DISCFIBERSCONTROL MATLAB code for ea_discfiberscontrol.fig
%      EA_DISCFIBERSCONTROL, by itself, creates a new EA_DISCFIBERSCONTROL or raises the existing
%      singleton*.
%
%      H = EA_DISCFIBERSCONTROL returns the handle to a new EA_DISCFIBERSCONTROL or the handle to
%      the existing singleton*.
%
%      EA_DISCFIBERSCONTROL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in EA_DISCFIBERSCONTROL.M with the given input arguments.
%
%      EA_DISCFIBERSCONTROL('Property','Value',...) creates a new EA_DISCFIBERSCONTROL or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ea_discfiberscontrol_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ea_discfiberscontrol_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ea_discfiberscontrol

% Last Modified by GUIDE v2.5 30-Jan-2019 11:26:57

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ea_discfiberscontrol_OpeningFcn, ...
                   'gui_OutputFcn',  @ea_discfiberscontrol_OutputFcn, ...
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


% --- Executes just before ea_discfiberscontrol is made visible.
function ea_discfiberscontrol_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ea_discfiberscontrol (see VARARGIN)

% Choose default command line output for ea_discfiberscontrol
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ea_discfiberscontrol wait for user response (see UIRESUME)
% uiwait(handles.discfiberscontrol);

resultfig = varargin{1};
showfibersset = getappdata(resultfig, 'showfibersset');
pospredthreshold = round(getappdata(resultfig, 'pospredthreshold')*100);
negpredthreshold = round(getappdata(resultfig, 'negpredthreshold')*100);

% Set current figure to draw slider
set(0, 'CurrentFigure', handles.discfiberscontrol);
set(handles.discfiberscontrol, 'Name', '');

% Bind showfibersset setting changed function
set(handles.showfiberssetpanel, 'SelectionChangedFcn', {@showfiberssetChanged, handles, resultfig});

%% threshold slider:
jSlider{1} = javax.swing.JSlider(0,100);
posthreshLabelPos = get(handles.posthresholdLabel, 'Position');
javacomponent(jSlider{1},[posthreshLabelPos(1)-6,posthreshLabelPos(2)-45,300,45]);
set(jSlider{1}, 'Name', 'posthresholdSlider', 'Value', pospredthreshold,...
    'Background', java.awt.Color(1,1,1),...
    'MajorTickSpacing', 0.1, 'PaintLabels',true);
set(handles.posthresholdValue,'String',sprintf('%d',pospredthreshold));
hjSlider{1} = handle(jSlider{1}, 'CallbackProperties');
set(hjSlider{1}, 'MouseReleasedCallback', {@sliderThresholdChange, handles, resultfig});
set(hjSlider{1}, 'StateChangedCallback', {@sliderThresholdChangeTxt, handles});

jSlider{2} = javax.swing.JSlider(0,100);
negthreshLabelPos = get(handles.negthresholdLabel, 'Position');
javacomponent(jSlider{2},[negthreshLabelPos(1)-6,negthreshLabelPos(2)-45,300,45]);
set(jSlider{2}, 'Name', 'negthresholdSlider', 'Value', negpredthreshold,...
    'Background', java.awt.Color(1,1,1),...
    'MajorTickSpacing', 0.1, 'PaintLabels',true);
set(handles.negthresholdValue,'String',sprintf('%d',negpredthreshold));
hjSlider{2} = handle(jSlider{2}, 'CallbackProperties');
set(hjSlider{2}, 'MouseReleasedCallback', {@sliderThresholdChange, handles, resultfig});
set(hjSlider{2}, 'StateChangedCallback', {@sliderThresholdChangeTxt, handles});

setappdata(handles.discfiberscontrol, 'hjSlider', hjSlider);

switch showfibersset
    case 'positive'
        set(handles.showposonly, 'Value', 1);
        hjSlider{1}.setEnabled(1);
        hjSlider{2}.setEnabled(0);
    case 'negative'
        set(handles.shownegonly, 'Value', 1);
        hjSlider{1}.setEnabled(0);
        hjSlider{2}.setEnabled(1);
    case 'both'
        set(handles.showboth, 'Value', 1);
        hjSlider{1}.setEnabled(1);
        hjSlider{2}.setEnabled(1);
end


% --- Outputs from this function are returned to the command line.
function varargout = ea_discfiberscontrol_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


function showfiberssetChanged(hObject, eventdata, handles, resultfig)
% Update showfibersset setting
hjSlider = getappdata(handles.discfiberscontrol, 'hjSlider');
switch get(get(handles.showfiberssetpanel, 'SelectedObject'), 'Tag')
    case 'showposonly'
        showfibersset = 'positive';
        hjSlider{1}.setEnabled(1);
        hjSlider{2}.setEnabled(0);
    case 'shownegonly'
        showfibersset = 'negative';
        hjSlider{1}.setEnabled(0);
        hjSlider{2}.setEnabled(1);
    case 'showboth'
        showfibersset = 'both';
        hjSlider{1}.setEnabled(1);
        hjSlider{2}.setEnabled(1);
end

% predthreshold values remain the same as in the resultfig
pospredthreshold = getappdata(resultfig, 'pospredthreshold');
negpredthreshold = getappdata(resultfig, 'negpredthreshold');

updateDiscFibers(resultfig, showfibersset, pospredthreshold, negpredthreshold)


function sliderThresholdChangeTxt(hObject, eventdata, handles)
slider = hObject;

if slider.Enabled
    switch slider.Name
        case 'posthresholdSlider'
            set(handles.posthresholdValue,'String',sprintf('%d',slider.Value));
        case 'negthresholdSlider'
             set(handles.negthresholdValue,'String',sprintf('%d',slider.Value));
    end
end


function sliderThresholdChange(hObject, eventdata, handles, resultfig)
slider = hObject;

if slider.Enabled
    sliderThresholdChangeTxt(hObject, eventdata, handles);
    
    % showfibersset setting from the control figure
    switch get(get(handles.showfiberssetpanel, 'SelectedObject'), 'Tag')
        case 'showposonly'
            showfibersset = 'positive';
        case 'shownegonly'
            showfibersset = 'negative';
        case 'showboth'
            showfibersset = 'both';
    end

    switch slider.Name
        case 'posthresholdSlider'
            % negpredthreshold value remains the same as in the resultfig
            pospredthreshold = slider.Value/100;
            negpredthreshold = getappdata(resultfig, 'negpredthreshold');
        case 'negthresholdSlider'
            % pospredthreshold value remains the same as in the resultfig
            pospredthreshold = getappdata(resultfig, 'pospredthreshold');
            negpredthreshold = slider.Value/100;
    end

    updateDiscFibers(resultfig, showfibersset, pospredthreshold, negpredthreshold);
    set(0, 'CurrentFigure', handles.discfiberscontrol);
end


function updateDiscFibers(resultfig, showfibersset, pospredthreshold, negpredthreshold)
% Delete previous discfibers and colorbar
delete(getappdata(resultfig, 'discfibers'));
delete(getappdata(resultfig, 'cbfig'));

% Set current figure to update discfibers
set(0, 'CurrentFigure', resultfig);

% vals and fibcell to be trimmed for visualization
tvals = getappdata(resultfig, 'vals');
tfibcell = getappdata(resultfig, 'fibcell');

% Determine threshold
posits = getappdata(resultfig, 'posits');
negits = getappdata(resultfig, 'negits');
posthresh=posits(round(length(posits)*pospredthreshold));
negthresh=negits(round(length(negits)*negpredthreshold));

switch showfibersset
    case 'positive'
        negthresh = negits(1)-eps;
        disp(['Fiber colors: Positive (T = ',num2str(posthresh),' ~ ',num2str(posits(1)), ')']);
    case 'negative'
        posthresh = posits(1)+eps;
        disp(['Fiber colors: Negative (T = ',num2str(negits(1)),' ~ ',num2str(negthresh), ')']);
    case 'both'
        disp(['Fiber colors: Positive (T = ',num2str(posthresh),' ~ ',num2str(posits(1)), ...
          '); Negative (T = ',num2str(negits(1)),' ~ ',num2str(negthresh),').']);
end

% Remove tvals and fibers outside the thresholding range
remove=logical(logical(tvals<posthresh) .* logical(tvals>negthresh));
tvals(remove)=[];
tfibcell(remove)=[];

% Rescale positive/negative tvals to [0 1]/[-1 0]
tvalsRescale = tvals;
tvalsRescale(tvals>0)=ea_minmax(tvals(tvals>0));
tvalsRescale(tvals<0)=-ea_minmax(-tvals(tvals<0));

% Contruct colormap
colormap gray
map=ea_redblue(1024);
fibcolorInd=tvalsRescale*(size(map,1)/2-0.5);
fibcolorInd=fibcolorInd+(size(map,1)/2+0.5);

% Set alphas of fibers with light color to 0
colorbarThreshold = 0.60; % Percentage of the pos/neg color to be kept
negUpperBound=ceil(size(map,1)/2*colorbarThreshold);
poslowerBound=floor((size(map,1)-size(map,1)/2*colorbarThreshold));
alphas=zeros(size(fibcolorInd,1),1);
switch showfibersset
    case 'positive'
        alphas(round(fibcolorInd)>=poslowerBound) = 1;
    case 'negative'
        alphas(round(fibcolorInd)<=negUpperBound) = 1;
    case 'both'
        alphas(round(fibcolorInd)>=poslowerBound) = 1;
        alphas(round(fibcolorInd)<=negUpperBound) = 1;
end

alphas(round(fibcolorInd)>=poslowerBound) = 1;
fibalpha=mat2cell(alphas,ones(size(fibcolorInd,1),1));

% Plot fibers
h=streamtube(tfibcell,0.2);
nones=repmat({'none'},size(fibcolorInd));
[h.EdgeColor]=nones{:};

% Calulate fiber colors
colors=map(round(fibcolorInd),:);
fibcolor=mat2cell(colors,ones(size(fibcolorInd)));

% Set fiber colors and alphas
[h.FaceColor]=fibcolor{:};
[h.FaceAlpha]=fibalpha{:};

setappdata(resultfig, 'discfibers', h);

% Set colorbar tick positions and labels
cbvals = tvals(logical(alphas));
% cbvals=tvalsRescale(logical(alphas));
switch showfibersset
    case 'positive'
        cbmap = map(ceil(length(map)/2+0.5):end,:);
        tick = [poslowerBound, length(map)] - floor(length(map)/2) ;
        poscbvals = sort(cbvals(cbvals>0));
        ticklabel = [poscbvals(1), poscbvals(end)];
        ticklabel = arrayfun(@(x) num2str(x,'%.2f'), ticklabel, 'Uni', 0);
    case 'negative'
        cbmap = map(1:floor(length(map)/2-0.5),:);
        tick = [1, negUpperBound];
        negcbvals = sort(cbvals(cbvals<0));
        ticklabel = [negcbvals(1), negcbvals(end)];
        ticklabel = arrayfun(@(x) num2str(x,'%.2f'), ticklabel, 'Uni', 0);
    case 'both'
        cbmap = map;
        tick = [1, negUpperBound, poslowerBound, length(map)];
        poscbvals = sort(cbvals(cbvals>0));
        negcbvals = sort(cbvals(cbvals<0));
        ticklabel = [min(cbvals), negcbvals(end), poscbvals(1), max(cbvals)];
        ticklabel = arrayfun(@(x) num2str(x,'%.2f'), ticklabel, 'Uni', 0);
end

% Plot colorbar
cbfig = ea_plot_colorbar(cbmap, [], 'v', '', tick, ticklabel);
set(cbfig, 'NumberTitle', 'off');
setappdata(resultfig, 'cbfig', cbfig);

% Update appdata in resultfig
setappdata(resultfig, 'showfibersset', showfibersset);
setappdata(resultfig, 'pospredthreshold', pospredthreshold);
setappdata(resultfig, 'negpredthreshold', negpredthreshold);
