function reply = buttondlg(headerStr,optionStr,defaultValue)
% function reply = buttondlg(headerStr,optionStr,<defaultValue>)
%
% Button Dialog Box
%
%  Allows the user to toggle and select options in
%  the cell array string 'optionStr'.  Reply is a
%  boolean vector with length(optionStr) that indicates
%  selected options. reply is empty if the user
%  chooses to cancel.
%
%  Example:
%  reply = buttondlg('pick it',{'this','that','the other'})
%
%  Example:
%  reply = buttondlg('pick it',{'this','that','the other'},[0 1 0])
%
%4/11/98 gmb   Wrote it
% 15/07/02  fwc     it will now make multiple columns of options if the number is large
% 2002.07.24  rfd & fwc: fixed minor bug where short checkbox labels were invisible.
% 2010/08/27 jb adapt widht of each column to max length string

OptionsPerColumn=25; % max number of options in one column

if nargin<2
    disp('Error: "bottondlg" requires two inputs');
    return
end
if nargin<3
  defaultValue = zeros(1,length(optionStr));
else
  if length(defaultValue) < length(optionStr)
    defaultValue(end+1:length(optionStr)) = 1;
  end
end

if isunix
    fontSize = 12;
else
    fontSize = 11;
end
fontName = 'Helvetica';
margin = 5;
minColumnWidth = 150;
minFigureWidth = 250;
checkBoxWidth = 20;

%open the figure
hFigure = figure('name',char(headerStr),'MenuBar','none','NumberTitle','off');
figurePosition = get(hFigure,'Position');
%get the MLR figure position
figloc = mrGetFigLoc('buttondlg');
if ~isempty(figloc)
  figurePosition(1:2) = figloc(1:2);
end  

if ischar(optionStr)
   optionStr = cellstr(optionStr);
end

nOptions = length(optionStr);

nMultiCols=ceil(nOptions/OptionsPerColumn);
nOptionsPerColumn=ceil(nOptions/nMultiCols);

maxOptionWidth = minColumnWidth*ones(1,nMultiCols);
%compute a width per column
for iColumn=1:nMultiCols
   %convert to char to find max length more easily
   tempOptionStrings = optionStr((iColumn-1)*nOptionsPerColumn+1:min(iColumn*nOptionsPerColumn,nOptions));
  for iOption = 1:length(tempOptionStrings)
      h = uicontrol(hFigure,'Style','text','String',tempOptionStrings{iOption},'FontSize',fontSize,'FontName',fontName);
      thisExtent = get(h,'extent');
      maxOptionWidth(iColumn) = max(maxOptionWidth(iColumn),thisExtent(3)+checkBoxWidth);
      delete(h);
  end
end
optionHeight = thisExtent(4);
maxOptionWidth = maxOptionWidth*max(1,minFigureWidth/sum(maxOptionWidth));

% set the figure position
figureWidth = sum(maxOptionWidth)+(nMultiCols+1)*margin+checkBoxWidth;
figureHeight = margin*2+(nOptionsPerColumn+1)*optionHeight;

figurePosition(3:4) = [figureWidth figureHeight];
%deal with multiple monitors
[monitorNumber,figurePosition] = getMonitorNumber(figurePosition,getMonitorPositions);
% monitorPositions = monitorPositions(monitorNumber,:);
% %make the figure large enough to get length of text boxes
% figurePosition(3) = monitorPositions(3);
% set(hFigure,'position',figurePosition);

set(hFigure,'Position',figurePosition);

%Display the checkboxes
backGroundColor = get(hFigure,'Color');
y0 = figurePosition(4) - margin - optionHeight;
cMultiCol=0;
hButton = zeros(1,nOptions);
for optionNum=1:nOptions
  if optionNum>(cMultiCol)*nOptionsPerColumn
      cMultiCol=cMultiCol+1;
      y=y0;
  end
  hButton(optionNum) = ...
      uicontrol('Style','checkbox',...
      'Value',defaultValue(optionNum),...
      'Units','normalized',...
      'String',optionStr{optionNum},...
      'BackgroundColor',backGroundColor,...
      'Position',[(sum(maxOptionWidth(1:cMultiCol-1))+margin)/figureWidth  y/figureHeight...
                  maxOptionWidth(cMultiCol)/figureWidth  optionHeight/figureHeight],...
      'HorizontalAlignment','left',...
      'FontSize',fontSize);
  y = y-optionHeight;
end
if rem(nOptions,nOptionsPerColumn)
  y = y-optionHeight;
end
%Display the OK/Cancel buttons
buttonWidth = min(130,(figurePosition(3)-2*margin)/3);
intervalBetweenButtons = (figurePosition(3) - buttonWidth*2)/3;

uicontrol('Style','pushbutton',...
    'String','Cancel',...
    'units','normalized',...
    'Position',[intervalBetweenButtons/figureWidth y/figureHeight buttonWidth/figureWidth optionHeight/figureHeight],...
    'CallBack','uiresume',...
    'FontSize',fontSize,...
    'UserData','Cancel');

uicontrol('Style','pushbutton',...
    'String','OK',...
    'units','normalized',...
    'Position',[(intervalBetweenButtons*2+buttonWidth)/figureWidth y/figureHeight buttonWidth/figureWidth optionHeight/figureHeight],...
    'CallBack','uiresume',...
    'FontSize',fontSize,...
    'UserData','OK');

%let the user select some radio buttons and
%wait for a 'uiresume' callback from OK/Cancel
uiwait

%determine which button was hit.
response = get(gco,'UserData');

%gather the status of the radio buttons if 'OK' was
%selected.  Otherwise return empty matrix.
if strcmp(response,'OK')
  for optionNum=1:nOptions
      reply(optionNum)=get(hButton(optionNum),'Value');
  end
elseif strcmp(response,'Cancel') %here we differentiate between a Cancel press
    reply = [];
else                             %and the user closing the window with the left (right) top button
    reply = zeros(0,1);
end

if ishandle(hFigure)
  figurePosition = get(hFigure,'Position');
  mrSetFigLoc('buttondlg',figurePosition);
  close(hFigure)
  drawnow;
end
