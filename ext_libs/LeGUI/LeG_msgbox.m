function varargout=LeG_msgbox(varargin)
%MSGBOX Message box.
%   msgbox(Message) creates a message box that automatically wraps
%   Message to fit an appropriately sized Figure.  Message is a string
%   vector, string matrix or cell array.
%
%   msgbox(Message,Title) specifies the title of the message box.
%
%   msgbox(Message,Title,Icon) specifies which Icon to display in
%   the message box.  Icon is 'none', 'error', 'help', 'warn', or
%   'custom'. The default is 'none'.
%
%   msgbox(Message,Title,'custom',IconData,IconCMap) defines a
%   customized icon.  IconData contains image data defining the icon;
%   IconCMap is the colormap used for the image.
%
%   msgbox(Message,...,CreateMode) specifies whether a message box is modal
%   or non-modal. Valid values for CreateMode are 'modal', 'non-modal', and
%   'replace'. If CreateMode is 'modal' or 'replace', the first available
%   message box with the specified title is updated to reflect the new
%   properties of the message box. All other such message boxes are deleted.
%   If CreateMode is 'non-modal', the message-box is not replaced and a new
%   handle is created. The default value for CreateMode is 'non-modal'.
%
%   CreateMode may also be a structure with fields WindowStyle and
%   Interpreter.  WindowStyle may be any of the values above.
%   Interpreter may be 'tex' or 'none'.  The default value for the
%   Interpreter is 'none';
%
%   h = msgbox(...) returns the handle of the box in h.
%
%   To make msgbox block execution until the user responds, include the
%   string 'modal' in the input argument list and wrap the call to
%   msgbox with UIWAIT.
%
%   Examples:
%       %An example which blocks execution until the user responds:
%       uiwait(msgbox('String','Title','modal'));
%
%       %An example using a custom Icon is:
%       Data=1:64;Data=(Data'*Data)/64;
%       h=msgbox('String','Title','custom',Data,hot(64))
%
%       %An example which reuses the existing msgbox window:
%       CreateStruct.WindowStyle='replace';
%       CreateStruct.Interpreter='tex';
%       h=msgbox('X^2 + Y^2','Title','custom',Data,hot(64),CreateStruct);
%
%   See also DIALOG, ERRORDLG, HELPDLG, INPUTDLG, LISTDLG,
%            QUESTDLG, TEXTWRAP, UIWAIT, WARNDLG.

%  Copyright 1984-2010 The MathWorks, Inc.


%%%%%%%%%%%%%%%%%%%%
%%% Nargin Check %%%
%%%%%%%%%%%%%%%%%%%%
narginchk(1,6);
nargoutchk(0,1);

if nargin > 0
    [varargin{:}] = convertStringsToChars(varargin{:});
end

inputStr = varargin{1};
BodyTextString = LeG_dialogCellstrHelper(inputStr);

% setup defaults
TitleString=' ';
IconString ='none';
IconData   =[];
IconCMap   =[];


createArg = '';
if nargin > 1
    createArg = varargin{nargin};
end

[Flag,CreateMode,Interpreter]=InternalCreateFlag(createArg);

% Do a check on the Interpreter property upfront.
if (any(strcmpi(Interpreter, {'latex', 'tex', 'none'})) ~= 1)
    error(message('MATLAB:msgbox:interpreter'));
end


switch nargin
    case 2
        if ~Flag
            TitleString=varargin{2};
        end
    case 3
        TitleString=varargin{2};
        if ~Flag
            IconString=varargin{3};
        end
    case 4
        TitleString=varargin{2};
        IconString=varargin{3};
        if ~Flag
            IconData = varargin{4};
        end
    case 5
        if Flag
            error(message('MATLAB:msgbox:colormap'));
        end
        TitleString=varargin{2};
        IconString=varargin{3};
        if ~strcmpi(IconString,'custom')
            warning(message('MATLAB:msgbox:customicon'));
            IconString='custom';
        end
        IconData=varargin{4};
        IconCMap=varargin{5};
    case 6
        TitleString=varargin{2};
        IconString=varargin{3};
        IconData=varargin{4};
        IconCMap=varargin{5};
end

IconString=lower(IconString);
switch(IconString)
    case {'custom'}
        % check for icon data
        if isempty(IconData)
            error(message('MATLAB:msgbox:icondata'))
        end
        if ~isnumeric(IconData)
            error(message('MATLAB:msgbox:IncorrectIconDataType'))
        end
        if ~isnumeric(IconCMap)
            error(message('MATLAB:msgbox:IncorrectIconColormap'))
        end
    case {'none','help','warn','error'}
        % icon String OK
    otherwise
        warning(message('MATLAB:msgbox:iconstring'));
        IconString='none';
end

Black = [0 0 0];

%%%%%%%%%%%%%%%%%%%%%
%%% Set Positions %%%
%%%%%%%%%%%%%%%%%%%%%
DefFigPos=get(0,'DefaultFigurePosition');

MsgOff=7;
IconWidth = 32 * 72/get(groot,'ScreenPixelsPerInch');
IconHeight = 32 * 72/get(groot,'ScreenPixelsPerInch');

if strcmp(IconString,'none')
    FigWidth=125;
    if(~isunix)
        % Figure width for windows
        FigWidth=150;
    end
    MsgTxtWidth=FigWidth-2*MsgOff;
else
    FigWidth=190;
    MsgTxtWidth=FigWidth-2*MsgOff-IconWidth;
end
FigHeight=50;
DefFigPos(3:4)=[FigWidth FigHeight];

OKWidth=40;
OKHeight=17;
OKXOffset=(FigWidth-OKWidth)/2;
OKYOffset=MsgOff;


MsgTxtXOffset=MsgOff;
MsgTxtYOffset=MsgOff+OKYOffset+OKHeight;
MsgTxtHeight=FigHeight-MsgOff-MsgTxtYOffset;
MsgTxtForeClr=Black;

IconXOffset=MsgTxtXOffset;
IconYOffset=FigHeight-MsgOff-IconHeight;

%%%%%%%%%%%%%%%%%%%%%
%%% Create MsgBox %%%
%%%%%%%%%%%%%%%%%%%%%

figureHandle=[];

% See if a modal or replace dialog already exists and delete all of its
% children
MsgboxTag = ['Msgbox_', TitleString];
if ~strcmp(CreateMode,'non-modal')
    TempHide=get(0,'ShowHiddenHandles');
    set(0,'ShowHiddenHandles','on');
    OldFig=findobj(0,'Type','figure','Tag',MsgboxTag,'Name',TitleString);
    set(0,'ShowHiddenHandles',TempHide);
    if ~isempty(OldFig)
        figureHandle=OldFig;
        if length(OldFig)>1
            figureHandle=OldFig(1);
            close(OldFig(2:end));
            OldFig(2:end)=[];  %#ok
        end % if length
        CurPos=get(figureHandle,'Position');
        CurPos(3:4)=[FigWidth FigHeight];
        set(figureHandle,'Position',CurPos);
        BoxChildren=get(figureHandle,'Children');
        delete(BoxChildren);
        figure(figureHandle);
    end
end

if strcmpi(CreateMode,'modal')
    WindowStyle='modal';
else
    WindowStyle='normal';
end

if isempty(figureHandle)
    figureHandle=dialog(                                ...
        'Name'            ,TitleString             , ...
        'Pointer'         ,'arrow'                 , ...
        'Units'           ,'points'                , ...
        'Visible'         ,'off'                   , ...
        'KeyPressFcn'     ,@doKeyPress             , ...
        'WindowStyle'     ,WindowStyle             , ...
        'Toolbar'         ,'none'                  , ...
        'HandleVisibility','on'                    , ...
        'Tag'             ,MsgboxTag                 ...
        );
    % should this be 'on' to match the case below?
    %'HandleVisibility','callback'              , ...

else
    set(figureHandle,   ...
        'WindowStyle'     ,WindowStyle, ...
        'HandleVisibility','on'         ...
        );
end 

FigColor=get(figureHandle,'Color');

MsgTxtBackClr=FigColor;

Font.FontUnits='points';
Font.FontSize=get(0,'FactoryUicontrolFontSize');
Font.FontName=get(0,'FactoryUicontrolFontName');
Font.FontWeight=get(figureHandle,'DefaultUicontrolFontWeight');

StFont = Font;
StFont.FontSize = 12;
StFont.FontWeight=get(figureHandle, 'DefaultTextFontWeight');

okPos = [ OKXOffset OKYOffset OKWidth OKHeight ];
if isempty(TitleString)
    OKHandle=uicontrol(figureHandle                             , ...
        Font                                                    , ...
        'Style'              ,'pushbutton'                      , ...
        'Units'              ,'points'                          , ...
        'Position'           , okPos                            , ...
        'Callback'           ,'delete(gcbf)'                    , ...
        'KeyPressFcn'        ,@doKeyPress                       , ...
        'String'             ,getString(message('MATLAB:uistring:popupdialogs:OK'))                              , ...
        'HorizontalAlignment','center'                          , ...
        'Tag'                ,'OKButton'                          ...
        );
else
    OKHandle=uicontrol(figureHandle                             , ...
        Font                                                    , ...
        'Style'              ,'pushbutton'                      , ...
        'Units'              ,'points'                          , ...
        'Position'           , okPos                            , ...
        'Callback'           ,{@openPNG,TitleString}            , ...
        'KeyPressFcn'        ,@doKeyPress                       , ...
        'String'             ,'Example'                         , ...
        'HorizontalAlignment','center'                          , ...
        'Tag'                ,'OKButton'                          ...
        );
end

msgPos = [ MsgTxtXOffset MsgTxtYOffset MsgTxtWidth MsgTxtHeight ];
MsgHandle=uicontrol(figureHandle         , ...
    StFont                               , ...
    'Style'              ,'text'         , ...
    'Units'              ,'points'       , ...
    'Position'           , msgPos        , ...
    'String'             ,' '            , ...
    'Tag'                ,'MessageBox'   , ...
    'HorizontalAlignment','left'         , ...
    'BackgroundColor'    ,MsgTxtBackClr  , ...
    'ForegroundColor'    ,MsgTxtForeClr    ...
    );

[WrapString,NewMsgTxtPos]=textwrap(MsgHandle,BodyTextString,MsgTxtWidth);
delete(MsgHandle);

% place an axes for the messge string (use an axes so we can get
% latex interpreter if required
AxesHandle=axes( ...
    'Parent'             ,figureHandle  , ...
    'Position'           ,[0 0 1 1]     , ...
    'Visible'            ,'off'           ...
    );

texthandle = text( ...
    'Parent'              ,AxesHandle                        , ...
    'Units'               ,'points'                          , ...
    'String'              ,WrapString                        , ...
    'Color'               ,get(OKHandle,'ForegroundColor')   , ...
    StFont                                                   , ...
    'HorizontalAlignment' ,'left'                            , ...
    'VerticalAlignment'   ,'bottom'                          , ...
    'Interpreter'         ,Interpreter                       , ...
    'Tag'                 ,'MessageBox'                        ...
    );

textExtent = get(texthandle, 'Extent');

%textExtent and extent from uicontrol are not the same. For window, extent from uicontrol is larger
%than textExtent. But on Mac, it is reverse. Pick the max value.
MsgTxtWidth=max([MsgTxtWidth NewMsgTxtPos(3) textExtent(3)]);
MsgTxtHeight=max([MsgTxtHeight NewMsgTxtPos(4) textExtent(4)]);

if ~strcmp(IconString,'none')
    MsgTxtXOffset=IconXOffset+IconWidth+MsgOff;
    FigWidth=MsgTxtXOffset+MsgTxtWidth+MsgOff;
    % Center Vertically around icon
    if IconHeight>MsgTxtHeight
        IconYOffset=OKYOffset+OKHeight+MsgOff;
        MsgTxtYOffset=IconYOffset+(IconHeight-MsgTxtHeight)/2;
        FigHeight=IconYOffset+IconHeight+MsgOff;
        % center around text
    else
        MsgTxtYOffset=OKYOffset+OKHeight+MsgOff;
        IconYOffset=MsgTxtYOffset+(MsgTxtHeight-IconHeight)/2;
        FigHeight=MsgTxtYOffset+MsgTxtHeight+MsgOff;
    end

else
    FigWidth=MsgTxtWidth+2*MsgOff;
    MsgTxtYOffset=OKYOffset+OKHeight+MsgOff;
    FigHeight=MsgTxtYOffset+MsgTxtHeight+MsgOff;
end % if ~strcmp

OKXOffset=(FigWidth-OKWidth)/2;
DefFigPos(3:4)=[FigWidth FigHeight];
DefFigPos = LeG_getnicedialoglocation(DefFigPos, get(figureHandle,'Units'));

% if there is a figure out there and it's modal, we need to be modal too
if ~isempty(gcbf) && strcmp(get(gcbf,'WindowStyle'),'modal')
    set(figureHandle,'WindowStyle','modal');
end

set(figureHandle,'Position',DefFigPos);

set(OKHandle,'Position',[OKXOffset OKYOffset OKWidth OKHeight]);

% calculate location for shadow box and put behind the button
LeG_setdefaultbutton(figureHandle, OKHandle);

txtPos = [ MsgTxtXOffset MsgTxtYOffset 0 ];
set(texthandle, 'Position'            ,txtPos);

if ~strcmp(IconString,'none')
    % create an axes for the icon
    iconPos = [IconXOffset IconYOffset IconWidth IconHeight];
    IconAxes=axes(                                   ...
        'Parent'          ,figureHandle               , ...
        'Units'           ,'points'                , ...
        'Position'        ,iconPos                 , ...
        'Tag'             ,'IconAxes'                ...
        );

    if ~strcmp(IconString,'custom')
        % Cases where IconString will be one of 'help','warn' or 'error'
        Img = setupStandardIcon(IconAxes, IconString);        
    else
        % place the icon - if this fails, rethrow the error
        % after deleting the figure
        try
            Img=image('CData',IconData,'Parent',IconAxes);
            set(figureHandle, 'Colormap', IconCMap);
        catch ex
            delete(figureHandle);
            rethrow(ex);
        end
    end
    if ~isempty(get(Img,'XData')) && ~isempty(get(Img,'YData'))
        set(IconAxes          , ...
            'XLim'            ,get(Img,'XData')+[-0.5 0.5], ...
            'YLim'            ,get(Img,'YData')+[-0.5 0.5]  ...
            );
    end

    set(IconAxes          , ...
        'Visible'         ,'off'       , ...
        'YDir'            ,'reverse'     ...
        );

end % if ~strcmp

% make sure we are on screen
movegui(figureHandle)

set(figureHandle,'HandleVisibility','callback','Visible','on');

% make sure the window gets drawn even if we are in a pause
drawnow

if nargout==1
    varargout{1}=figureHandle;
end

end

%%%%% InternalCreateFlag
function [Flag,CreateMode,Interpreter]=InternalCreateFlag(mode)
    Flag=0;
    CreateMode='non-modal';
    Interpreter='none';

    if isempty(mode)
        return;
    end

    if iscell(mode)
        mode=mode{:};
    end

    if isstruct(mode)

        if ~isfield(mode,'Interpreter') || ~isfield(mode,'WindowStyle')
            error(message('MATLAB:msgbox:InvalidInput'));
        end
        
        Interpreter=mode.Interpreter;
        mode=mode.WindowStyle;
    end

    if ~ischar(mode)
        return;
    end

    mode=lower(mode);
    switch(mode)
        case {'non-modal','modal','replace'}
         CreateMode = mode;
         Flag=1;
    end
end

%%%%% doKeyPress
function doKeyPress(obj, evd)
    switch(evd.Key)
        case {'return','space','escape'}
            delete(ancestor(obj,'figure'));
    end
end

function Img = setupStandardIcon(ax, iconName)
[iconData, alphaData] = matlab.ui.internal.dialog.DialogUtils.imreadDefaultIcon(iconName);  
Img=image('CData',iconData,'Parent',ax);
if ~isempty(alphaData)
    set(Img, 'AlphaData', alphaData)
end
end

function openPNG(varargin)
fH = figure('menubar','none','name','Example');
aH = axes('parent',fH);
if isdeployed
    switch varargin{3}
        case 'ACPCGUI'
            imshow('ACPCGUIExample.png','border','tight','interpolation','bilinear','parent',aH);
        case 'AssignElectrodesGUI'
            imshow('AssignElectrodesExample.png','border','tight','interpolation','bilinear','parent',aH);
    end
else
    switch varargin{3}
        case 'ACPCGUI'
            imshow('./icons/ACPCGUIExample.png','border','tight','interpolation','bilinear','parent',aH);
        case 'AssignElectrodesGUI'
            imshow('./icons/AssignElectrodesExample.png','border','tight','interpolation','bilinear','parent',aH);
    end
end
end
