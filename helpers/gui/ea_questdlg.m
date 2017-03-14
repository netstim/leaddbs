function ButtonName=ea_questdlg(Question,Title,Btn1,Btn2,Btn3,Default)
%QUESTDLG Question dialog box.
%  ButtonName = QUESTDLG(Question) creates a modal dialog box that
%  automatically wraps the cell array or string (vector or matrix)
%  Question to fit an appropriately sized window.  The name of the
%  button that is pressed is returned in ButtonName.  The Title of
%  the figure may be specified by adding a second string argument:
%
%    ButtonName = questdlg(Question, Title)
%
%  Question will be interpreted as a normal string.
%
%  QUESTDLG uses UIWAIT to suspend execution until the user responds.
%
%  The default set of buttons names for QUESTDLG are 'Yes','No' and 'Cancel'.
%  The default answer for the above calling syntax is 'Yes'.
%  This can be changed by adding a third argument which specifies the
%  default Button:
%
%    ButtonName = questdlg(Question, Title, 'No')
%
%  Up to 3 custom button names may be specified by entering
%  the button string name(s) as additional arguments to the function
%  call.  If custom button names are entered, the default button
%  must be specified by adding an extra argument, DEFAULT, and
%  setting DEFAULT to the same string name as the button you want
%  to use as the default button:
%
%    ButtonName = questdlg(Question, Title, Btn1, Btn2, DEFAULT);
%
%  where DEFAULT is set to Btn1.  This makes Btn1 the default answer.
%  If the DEFAULT string does not match any of the button string names,
%  a warning message is displayed.
%
%  To use TeX interpretation for the Question string, a data
%  structure must be used for the last argument, i.e.
%
%    ButtonName = questdlg(Question, Title, Btn1, Btn2, OPTIONS);
%
%  The OPTIONS structure must include the fields Default and Interpreter.
%  Interpreter may be 'none' or 'tex' and Default is the default button
%  name to be used.
%
%  If the dialog is closed without a valid selection, the return value
%  is empty.
%
%  Example:
%
%  ButtonName = questdlg('What is your favorite color?', ...
%                        'Color Question', ...
%                        'Red', 'Green', 'Blue', 'Green');
%  switch ButtonName,
%    case 'Red',
%     disp('Your favorite color is Red');
%    case 'Blue',
%     disp('Your favorite color is Blue.')
%     case 'Green',
%      disp('Your favorite color is Green.');
%  end % switch
%
%  See also DIALOG, ERRORDLG, HELPDLG, INPUTDLG, LISTDLG,
%    MSGBOX, WARNDLG, FIGURE, TEXTWRAP, UIWAIT, UIRESUME.


%  Copyright 1984-2014 The MathWorks, Inc.



if nargin<1
    error(message('MATLAB:questdlg:TooFewArguments'));
end

Interpreter='none';
% Question = dialogCellstrHelper(Question);
needsLookup = false;

%%%%%%%%%%%%%%%%%%%%%
%%% General Info. %%%
%%%%%%%%%%%%%%%%%%%%%
Black      =[0       0        0      ]/255;
% LightGray  =[192     192      192    ]/255;
% LightGray2 =[160     160      164    ]/255;
% MediumGray =[128     128      128    ]/255;
% White      =[255     255      255    ]/255;

%%%%%%%%%%%%%%%%%%%%
%%% Nargin Check %%%
%%%%%%%%%%%%%%%%%%%%
if nargout > 1
    error(message('MATLAB:questdlg:WrongNumberOutputs'));
end
if nargin == 1
    Title = ' ';
end
if nargin <= 2,
    Default = 'Yes'; 
    needsLookup = true;
end
if nargin == 3,
    Default = Btn1;
end
if nargin <= 3,
    Btn1 = 'Yes'; 
    Btn2 = 'No';
    Btn3 = 'Cancel';
    NumButtons = 3;
    needsLookup = true;
end
if nargin == 4,
    Default=Btn2;
    Btn2 = [];
    Btn3 = [];
    NumButtons = 1;
end
if nargin == 5
    Default = Btn3;
    Btn3 = [];
    NumButtons = 2;
end
if nargin == 6
    NumButtons = 3;
end
if nargin > 6
    error(message('MATLAB:questdlg:TooManyInputs'));
    NumButtons = 3;
end

if isstruct(Default),
    Interpreter = Default.Interpreter;
    Default = Default.Default;
end


%%%%%%%%%%%%%%%%%%%%%%%
%%% Create QuestFig %%%
%%%%%%%%%%%%%%%%%%%%%%%
FigPos    = get(0,'DefaultFigurePosition');
FigPos(3) = 267;
FigPos(4) =  70;
FigPos    = ea_getnicedialoglocation(FigPos, get(0,'DefaultFigureUnits'));

QuestFig=dialog(                                    ...
    'Visible'         ,'off'                      , ...
    'Name'            ,Title                      , ...
    'Pointer'         ,'arrow'                    , ...
    'Position'        ,FigPos                     , ...
    'KeyPressFcn'     ,@doFigureKeyPress          , ...
    'IntegerHandle'   ,'off'                      , ...
    'WindowStyle'     ,'normal'                   , ...
    'HandleVisibility','callback'                 , ...
    'CloseRequestFcn' ,@doDelete                  , ...
    'Tag'             ,Title                        ...
    );


%%%%%%%%%%%%%%%%%%%%%
%%% Set Positions %%%
%%%%%%%%%%%%%%%%%%%%%
DefOffset  =10;

IconWidth = 32;
IconHeight = 32;
IconXOffset=DefOffset;
IconYOffset=FigPos(4)-DefOffset-IconHeight;  %#ok
IconCMap=[Black;get(QuestFig,'Color')];  %#ok

DefBtnWidth =56;
BtnHeight   =22;

BtnYOffset=DefOffset;

BtnWidth=DefBtnWidth;

ExtControl=uicontrol(QuestFig   , ...
    'Style'    ,'pushbutton', ...
    'String'   ,' '          ...
    );

btnMargin=1.4;
set(ExtControl,'String',Btn1);
BtnExtent=get(ExtControl,'Extent');
BtnWidth=max(BtnWidth,BtnExtent(3)+8);
if NumButtons > 1
    set(ExtControl,'String',Btn2);
    BtnExtent=get(ExtControl,'Extent');
    BtnWidth=max(BtnWidth,BtnExtent(3)+8);
    if NumButtons > 2
        set(ExtControl,'String',Btn3);
        BtnExtent=get(ExtControl,'Extent');
        BtnWidth=max(BtnWidth,BtnExtent(3)*btnMargin);
    end
end
BtnHeight = max(BtnHeight,BtnExtent(4)*btnMargin);

delete(ExtControl);

MsgTxtXOffset=IconXOffset+IconWidth;

FigPos(3)=max(FigPos(3),MsgTxtXOffset+NumButtons*(BtnWidth+2*DefOffset));
set(QuestFig,'Position',FigPos);

BtnXOffset=zeros(NumButtons,1);

if NumButtons==1,
    BtnXOffset=(FigPos(3)-BtnWidth)/2;
elseif NumButtons==2,
    BtnXOffset=[MsgTxtXOffset
        FigPos(3)-DefOffset-BtnWidth];
elseif NumButtons==3,
    BtnXOffset=[MsgTxtXOffset
        0
        FigPos(3)-DefOffset-BtnWidth];
    BtnXOffset(2)=(BtnXOffset(1)+BtnXOffset(3))/2;
end

MsgTxtYOffset=DefOffset+BtnYOffset+BtnHeight;
% Calculate current msg text width and height. If negative,
% clamp it to 1 since its going to be recalculated/corrected later
% based on the actual msg string
MsgTxtWidth=max(1, FigPos(3)-DefOffset-MsgTxtXOffset-IconWidth);
MsgTxtHeight=max(1, FigPos(4)-DefOffset-MsgTxtYOffset);

MsgTxtForeClr=Black;
MsgTxtBackClr=get(QuestFig,'Color');

CBString='uiresume(gcbf)';
DefaultValid = false;
DefaultWasPressed = false;
BtnHandle = cell(NumButtons, 1);
DefaultButton = 0;

% Check to see if the Default string passed does match one of the
% strings on the buttons in the dialog. If not, throw a warning.
for i = 1:NumButtons
    switch i
        case 1
            ButtonString=Btn1;
            ButtonTag='Btn1';
            if strcmp(ButtonString, Default)
                DefaultValid = true;
                DefaultButton = 1;
            end
            
        case 2
            ButtonString=Btn2;
            ButtonTag='Btn2';
            if strcmp(ButtonString, Default)
                DefaultValid = true;
                DefaultButton = 2;
            end
        case 3
            ButtonString=Btn3;
            ButtonTag='Btn3';
            if strcmp(ButtonString, Default)
                DefaultValid = true;
                DefaultButton = 3;
            end
    end
    
    if (needsLookup)
        buttonDisplayString = getString(message(['MATLAB:uistring:popupdialogs:' ButtonString]));
    else
        buttonDisplayString = ButtonString;
    end
    
    BtnHandle{i}=uicontrol(QuestFig            , ...
        'Style'              ,'pushbutton', ...
        'Position'           ,[ BtnXOffset(1) BtnYOffset BtnWidth BtnHeight ]           , ...
        'KeyPressFcn'        ,@doControlKeyPress , ...
        'Callback'           ,CBString    , ...
        'String'             ,buttonDisplayString, ...
        'HorizontalAlignment','center'    , ...
        'Tag'                ,ButtonTag     ...
        );
    
    setappdata(BtnHandle{i},'QuestDlgReturnName',ButtonString);   
end

if ~DefaultValid
    warnstate = warning('backtrace','off');
    warning(message('MATLAB:questdlg:StringMismatch'));
    warning(warnstate);
end

MsgHandle=uicontrol(QuestFig            , ...
    'Style'              ,'text'         , ...
    'Position'           ,[MsgTxtXOffset MsgTxtYOffset 0.95*MsgTxtWidth MsgTxtHeight ]              , ...
    'String'             ,{' '}          , ...
    'Tag'                ,'Question'     , ...
    'HorizontalAlignment','left'         , ...
    'FontWeight'         ,'bold'         , ...
    'BackgroundColor'    ,MsgTxtBackClr  , ...
    'ForegroundColor'    ,MsgTxtForeClr    ...
    );

[WrapString,NewMsgTxtPos]=textwrap(MsgHandle,Question,75);

% NumLines=size(WrapString,1);

AxesHandle=axes('Parent',QuestFig,'Position',[0 0 1 1],'Visible','off');

texthandle=text( ...
    'Parent'              ,AxesHandle                      , ...
    'Units'               ,'pixels'                        , ...
    'Color'               ,get(BtnHandle{1},'ForegroundColor')   , ...
    'HorizontalAlignment' ,'left'                          , ...
    'FontName'            ,get(BtnHandle{1},'FontName')    , ...
    'FontSize'            ,get(BtnHandle{1},'FontSize')    , ...
    'VerticalAlignment'   ,'bottom'                        , ...
    'String'              ,WrapString                      , ...
    'Interpreter'         ,Interpreter                     , ...
    'Tag'                 ,'Question'                        ...
    );

textExtent = get(texthandle, 'Extent');

% (g357851)textExtent and extent from uicontrol are not the same. For window, extent from uicontrol is larger
%than textExtent. But on Mac, it is reverse. Pick the max value.
MsgTxtWidth=max([MsgTxtWidth NewMsgTxtPos(3)+2 textExtent(3)]);
MsgTxtHeight=max([MsgTxtHeight NewMsgTxtPos(4)+2 textExtent(4)]);

MsgTxtXOffset=IconXOffset+IconWidth+DefOffset;
FigPos(3)=max(NumButtons*(BtnWidth+DefOffset)+DefOffset, ...
    MsgTxtXOffset+MsgTxtWidth+DefOffset);


% Center Vertically around icon
if IconHeight>MsgTxtHeight,
    IconYOffset=BtnYOffset+BtnHeight+DefOffset;
    MsgTxtYOffset=IconYOffset+(IconHeight-MsgTxtHeight)/2;
    FigPos(4)=IconYOffset+IconHeight+DefOffset;
    % center around text
else
    MsgTxtYOffset=BtnYOffset+BtnHeight+DefOffset;
    IconYOffset=MsgTxtYOffset+(MsgTxtHeight-IconHeight)/2;
    FigPos(4)=MsgTxtYOffset+MsgTxtHeight+DefOffset;
end

if NumButtons==1,
    BtnXOffset=(FigPos(3)-BtnWidth)/2;
elseif NumButtons==2,
    BtnXOffset=[(FigPos(3)-DefOffset)/2-BtnWidth
        (FigPos(3)+DefOffset)/2
        ];
    
elseif NumButtons==3,
    BtnXOffset(2)=(FigPos(3)-BtnWidth)/2;
    BtnXOffset=[BtnXOffset(2)-DefOffset-BtnWidth
        BtnXOffset(2)
        BtnXOffset(2)+BtnWidth+DefOffset
        ];
end

set(QuestFig ,'Position',ea_getnicedialoglocation(FigPos, get(QuestFig,'Units')));
assert(iscell(BtnHandle));
BtnPos=cellfun(@(bh)get(bh,'Position'), BtnHandle, 'UniformOutput', false);
BtnPos=cat(1,BtnPos{:});
BtnPos(:,1)=BtnXOffset;
BtnPos=num2cell(BtnPos,2);

assert(iscell(BtnPos));
cellfun(@(bh,pos)set(bh, 'Position', pos), BtnHandle, BtnPos, 'UniformOutput', false);

if DefaultValid
    ea_setdefaultbutton(QuestFig, BtnHandle{DefaultButton});
end

delete(MsgHandle);


set(texthandle, 'Position',[MsgTxtXOffset MsgTxtYOffset 0]);


IconAxes=axes(                                      ...
    'Parent'      ,QuestFig              , ...
    'Units'       ,'Pixels'              , ...
    'Position'    ,[IconXOffset IconYOffset IconWidth IconHeight], ...
    'NextPlot'    ,'replace'             , ...
    'Tag'         ,'IconAxes'              ...
    );

set(QuestFig ,'NextPlot','add');

[iconData, alphaData] = matlab.ui.internal.dialog.DialogUtils.imreadDefaultIcon('quest');  
Img=image('CData',iconData, 'AlphaData', alphaData, 'Parent',IconAxes);
set(IconAxes, ...
    'Visible','off'           , ...
    'YDir'   ,'reverse'       , ...
    'XLim'   ,get(Img,'XData')+[-0.5 0.5], ...
    'YLim'   ,get(Img,'YData')+[-0.5 0.5]  ...
    );

% make sure we are on screen
movegui(QuestFig)

set(QuestFig, 'Visible','on');
WinOnTop( QuestFig, true );
set(QuestFig ,'WindowStyle','modal');
drawnow;

if DefaultButton ~= 0
    uicontrol(BtnHandle{DefaultButton});
end

if ishghandle(QuestFig)
    % Go into uiwait if the figure handle is still valid.
    % This is mostly the case during regular use.
    c = matlab.ui.internal.dialog.DialogUtils.disableAllWindowsSafely();
    uiwait(QuestFig);
	delete(c);
end

% Check handle validity again since we may be out of uiwait because the
% figure was deleted.
if ishghandle(QuestFig)
    if DefaultWasPressed
        ButtonName=Default;
    else
        ButtonName = getappdata(get(QuestFig,'CurrentObject'),'QuestDlgReturnName');
    end
    doDelete;
else
    ButtonName='';
end
drawnow; % Update the view to remove the closed figure (g1031998)


    function doFigureKeyPress(obj, evd)  %#ok
    switch(evd.Key)
        case {'return','space'}
            if DefaultValid
                DefaultWasPressed = true;
                uiresume(gcbf);
            end
        case 'escape'
            doDelete
    end
    end

    function doControlKeyPress(obj, evd)  %#ok
    switch(evd.Key)
        case {'return'}
            if DefaultValid
                DefaultWasPressed = true;
                uiresume(gcbf);
            end
        case 'escape'
            doDelete
    end
    end

    function doDelete(varargin)
    delete(QuestFig);
    end
end
