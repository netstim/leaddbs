% mrParamsDialog.m
%
%      usage: mrParamsDialog(paramsInfo,<titleString>,<buttonWidth>,<callback>,<callbackArg>,<okCallback>,<cancelCallback>)
%  alt usage: You can also set variables explicitly, with the following syntax:
%             mrParamsDialog(paramsInfo,<titleString>,varargin)
%       e.g.: mrParamsDialog(paramsInfo,'This is the title','buttonWidth=1.5');
%             valid variable names are (buttonWidth,callback,callbackArg,okCallback,cancelCallback)
%             and also ignoreKeys (which keeps mrParamsDialog from allowing ESC to close it)
%         by: justin gardner, modified by julien besle
%       date: 03/13/07
%    purpose: creates a dialog for selection of parameters
%             see wiki for details
%        $Id$
%
function [params params2] = mrParamsDialog(varargin)

% check arguments
if nargin < 1
  help mrParamsDialog
  return
end

% if this is a cell array, it means to open up the figure
% using the variable name, default value pairs given
if iscell(varargin{1})
  % if empty paramsInfo just return
  if isempty(varargin{1}),params = [];return,end
  % otherwise init the dialog
  [params params2] = initFigure(varargin{1},varargin);
  % otherwise it is a callback
  
elseif isnumeric(varargin{1}) && isnumeric(varargin{2}) && isnumeric(varargin{3})% if it is 3 or 4 numbers then an entry field has been updated
  if length(varargin) == 3 
    buttonHandler(varargin{1},varargin{2},varargin{3});
  elseif length(varargin) == 4
    buttonHandler(varargin{1},varargin{2},varargin{3},varargin{4});
  end
  
else
  mrWarnDlg('(mrParamsDialog) unknown input parameter type');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% handle keyboard
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mrParamsKeyPressFcn(figHandle,keyEvent)

switch (keyEvent.Key)
  case {'return'}
   okHandler;
  case {'escape'}
   closeHandler;
  case {'f1'}
   helpHandler;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set up figure in first place
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [params params2] = initFigure(vars,otherParams)

params=[];
params2 = [];
% close any existing params window
closeHandler;

global gParams; %gParams will only have non-graphical info about the parameters, including handles to uicontrols and figures
% dParams will be a local structure with graphical information specific to the dialog box being drawn
% uiParams will be a local structure with general graphical information

% parse the input parameter string
[gParams.vars gParams.varinfo] = mrParamsParse(vars);

% parse the otherParams. The first otherParams is always the title
if length(otherParams) > 1
  titleStr = otherParams{2};
else
  titleStr = '';
end

% now if there is a second otherParams and it is a string, then
% it means that we have been passed in a "getArgs" type. Otherwise
% it is the old calling convention in which the order determines what
% variable is being set
buttonWidth = [];callback = [];callbackArg = [];okCallback = [];cancelCallback = [];modal=[];fullWidth = false;
if (length(otherParams) > 2)
  if ~ischar(otherParams{3})
    % get the arguments the old way, by order
    buttonWidth = otherParams{3};
    if length(otherParams) > 3,callback = otherParams{4}; end
    if length(otherParams) > 4,callbackArg = otherParams{5}; end
    if length(otherParams) > 5,okCallback = otherParams{6}; end
    if length(otherParams) > 6,cancelCallback = otherParams{7}; end
    if length(otherParams) > 7 &&ischar(otherParams{7})
      getArgs(otherParams(7:end));
    end
  else
    getArgs(otherParams(3:end));
  end
end

% default to using a modal dialog if there is no callback
if isempty(modal)
  if ~isempty(callback)
    modal = 0;
  else
    modal = 1;
  end
end

% get the figure
if ~isfield(gParams,'fignum') || (gParams.fignum == -1);
  % open figure, and turn off menu
  gParams.fignum = figure;
else
  figure(gParams.fignum);
end
gParams.figlocstr{1} = sprintf('mrParamsDialog_%s',fixBadChars(titleStr));
set(gParams.fignum,'MenuBar','none');
set(gParams.fignum,'NumberTitle','off');
set(gParams.fignum,'closeRequestFcn',@altCloseHandler);
if isempty(titleStr)
  set(gParams.fignum,'Name','Set parameters');
else
  set(gParams.fignum,'Name',titleStr);
end

% some basic info about location of controls
uiParams.maxFigHeightWidthRatio = 1.7; 
uiParams.minEntriesWidth = 200; %the minimum width of all the parameter entries
uiParams.maxEntriesWidth = 350; %the maximum width of all the parameter entries
uiParams.maxSingleFieldWidth = 100;
uiParams.minVarNameWidth = 70;
uiParams.margin = 3;
uiParams.fontsize = 12;
uiParams.fontname = 'Helvetica';
uiParams.leftMargin = 6;
uiParams.topMargin = 6;
uiParams.maxIncdecButtonWidth = 50;
uiParams.incdecMargin = 2;
%Matlab doesn't return the extent of multiline text, so we have to guess how much smaller the text height is 
%relative to the height of a textbox, in order to avoid making text boxes that are too large when text wraps
%(although there probably is a way to get this information)
uiParams.lineHeightRatio = .67; %approximate height ot multiline text. 

% button width is in fact a button scaling parameter    
if ~isempty(buttonWidth)
  uiParams.maxEntriesWidth = buttonWidth*uiParams.maxEntriesWidth;
  uiParams.maxSingleFieldWidth = buttonWidth*uiParams.maxSingleFieldWidth;
end


% Collect information for uicontrol
%initialize varname field info(first column)
uiParams.varNameWidth = 0;
uiParams.varName = repmat({''},1,length(gParams.vars));
%initialize entry field info (second column)
dParams.entryWidth = -inf(1,length(gParams.vars));
dParams.entryValue = zeros(1,length(gParams.vars));
dParams.entryString = repmat({{''}},1,length(gParams.vars));
dParams.testString = repmat({''},1,length(gParams.vars));
dParams.entryStyle = repmat({''},1,length(gParams.vars));
dParams.entryNumCols = ones(1,length(gParams.vars));
dParams.entryNumRows = ones(1,length(gParams.vars));
dParams.incdec = zeros(length(gParams.vars),2);
dParams.incdecType = repmat({'arrows'},1,length(gParams.vars));
dParams.numLines = ones(1,length(gParams.vars));
for i = 1:length(gParams.vars)
  if ~gParams.varinfo{i}.visible 
    dParams.numLines(i)=0; %no line for parameters that are not visible
  end
  %get variable name width
  uiParams.varName{i} = [gParams.vars{i}{1} '  ']; %add spaces on the right

  %get info about the entries
  switch(gParams.varinfo{i}.type)
    case 'pushbutton' 
      dParams.entryStyle{i} = 'pushbutton';
      if isfield(gParams.varinfo{i},'buttonString')
        dParams.entryString{i} = {['  ' gParams.varinfo{i}.buttonString '  ']};%we need to allow some space for the button features
      end
      dParams.testString(i) = dParams.entryString{i};

    case 'popupmenu' 
      dParams.entryStyle{i} = 'popupmenu';
      dParams.entryValue(i) = 1;
      dParams.entryString{i} = {gParams.varinfo{i}.value};
      %make up a string of Xs of lengh equal to the longest string in the menu list
      if strcmp(gParams.varinfo{i}.popuptype,'numeric')
	if iscell(dParams.entryString{i}{1})
	  thisEntryString = cell2mat(dParams.entryString{i}{1});
	else
	  thisEntryString = [dParams.entryString{i}{1}{:}];
	end
        charNum = length(num2str(max(thisEntryString)));
      else
        charNum = size(char(dParams.entryString{i}{1}),2);
      end 
      dParams.testString{i} =repmat('X',1,charNum+3);

    case 'statictext'
      dParams.entryString{i} = {gParams.varinfo{i}.value};
      dParams.testString(i) = dParams.entryString{i};
      dParams.entryStyle{i} = 'text';

    case 'checkbox'
      if isnumeric(gParams.varinfo{i}.value)
        dParams.entryValue(i) = gParams.varinfo{i}.value;
      else
	entryValue = str2num(gParams.varinfo{i}.value);
	if isempty(entryValue)
	  dParams.entryValue(i) = 0;
	else
	  dParams.entryValue(i) = entryValue;
	end
      end
      dParams.entryStyle{i} = 'checkbox';
      if isfield(gParams.varinfo{i},'buttonString')
        dParams.entryString{i} = {gParams.varinfo{i}.buttonString};%we need to allow some space for the button features
      end

    case 'string'
      dParams.entryString{i} = {gParams.varinfo{i}.value};
      if isfield(gParams.varinfo{i},'editable') && isequal(gParams.varinfo{i}.editable,0)
        dParams.entryStyle{i} = 'text';
        dParams.testString(i) = dParams.entryString{i};
      else
        dParams.entryStyle{i} = 'edit';
      end

    case 'numeric'
      dParams.entryString{i} = {thisNum2str(gParams.varinfo{i}.value)};
      if isfield(gParams.varinfo{i},'editable') && isequal(gParams.varinfo{i}.editable,0)
        dParams.entryStyle{i} = 'text';
        dParams.testString(i) = dParams.entryString{i};
      else
        dParams.entryStyle{i} = 'edit';
      end

    case 'array'
      dParams.entryString{i} = num2cell(gParams.varinfo{i}.value);
      if isfield(gParams.varinfo{i},'editable') && isequal(gParams.varinfo{i}.editable,0)
        dParams.entryStyle{i} = 'text';
      else
        dParams.entryStyle{i} = 'edit';
      end
      dParams.entryNumCols(i) = size(dParams.entryString{i},2);
      dParams.entryNumRows(i) = size(dParams.entryString{i},1);

    case 'stringarray'
      dParams.entryString{i} = gParams.varinfo{i}.value;
      if isfield(gParams.varinfo{i},'editable') && isequal(gParams.varinfo{i}.editable,0)
        dParams.entryStyle{i} = 'text';
      else
        dParams.entryStyle{i} = 'edit';
      end
      dParams.entryNumCols(i) = size(dParams.entryString{i},2);
      dParams.entryNumRows(i) = size(dParams.entryString{i},1);

    otherwise
       keyboard %unknown type...
%         dParams.entryString{i} = gParams.varinfo{i}.value;
%         dParams.entryStyle{i} = gParams.varinfo{i}.type;
  end


  if isfield(gParams.varinfo{i},'incdec')
    if (~isfield(gParams.varinfo{i},'enable') || isequal(gParams.varinfo{i}.enable,1))          % only display increment/decrement buttons if the control is editable and enabled
      dParams.incdec(i,:)=gParams.varinfo{i}.incdec;
      if isfield(gParams.varinfo{i},'incdecType')
        dParams.incdecType{i} = gParams.varinfo{i}.incdecType;
      end
    else
      gParams.varinfo{i}=rmfield(gParams.varinfo{i},'incdec');
    end
  end
end  

%optimize figure dimensions
[figurePosition,dParams,uiParams] = optimizeFigure(gParams.fignum,gParams.figlocstr{1},dParams,uiParams);
figWidth = figurePosition(3);
figHeight = figurePosition(4);

% make all entries full width if called for
if fullWidth
  dParams.entryWidth(:) = max(dParams.entryWidth);
end

%cap widths that are more than the max
dParams.entryWidth(dParams.entryWidth>dParams.allEntriesWidth)=dParams.allEntriesWidth;
for i = 1:length(dParams.entryStyle)
  %if it's gonna be an array compute the field width 
  if dParams.entryNumCols(i)>1 || dParams.entryNumRows(i)>1
      dParams.entryWidth(i) = min(dParams.allEntriesWidth/dParams.entryNumCols(i)-uiParams.margin,uiParams.maxSingleFieldWidth);
  %if it's not a popupmenu or a button, set the width to max
  elseif ~ismember(dParams.entryStyle,{'popupmenu','pushbutton'})
    dParams.entryWidth(i)=dParams.allEntriesWidth;
  end
end

%set control dimensions to normalized so that the figure resizes
set(gParams.fignum,'defaultUicontrolUnits','normalized'); 
%computing the normalized positions is taken care of by getUIControlPos


%--------------------------------- make entry buttons-----------------------------------
for i = 1:length(gParams.varinfo)
  % make ui for varname
  gParams.ui.varname(i) = makeUIcontrol(i,gParams.fignum,dParams,uiParams,'varname');
  % make ui for entry
  [gParams.ui.varentry{i} gParams.ui.incdec{i}{1} gParams.ui.incdec{i}{2}] =...
     makeUIcontrol(i,gParams.fignum,dParams,uiParams,'varentry');
  if isfield(gParams.varinfo{i},'incdec')
    switch(gParams.varinfo{i}.type)
      case 'string'
        values = mrStr2num(gParams.varinfo{i}.value);
      case {'numeric','array'}
        values = gParams.varinfo{i}.value;
    end
    for j=1:size(values,1)
      for k=1:size(values,2)
        enableArrows(values(j,k),i,j,k);
      end
    end
  end
end

% set ok and cancel callback
if ~isempty(okCallback)
  gParams.okCallback = okCallback;
  if ~isempty(callbackArg)
    gParams.callbackArg = callbackArg;
  end
end
if ~isempty(cancelCallback)
  gParams.cancelCallback = cancelCallback;
  if ~isempty(callbackArg)
    gParams.callbackArg = callbackArg;
  end
end

% for each value that controls another one, call the buttonHandler to
% set up the correct dependency
for i = 1:length(gParams.varinfo)
  if isfield(gParams.varinfo{i},'controls')
    buttonHandler(i,1,1);
  end
end

% check enable/visible options
for i = 1:length(gParams.varinfo)
  if isfield(gParams.varinfo{i},'enable') && isequal(gParams.varinfo{i}.enable,0)
      set(gParams.ui.varentry{i},'enable','off');
  end
  if isfield(gParams.varinfo{i},'visible') && isequal(gParams.varinfo{i}.visible,0)
    set(gParams.ui.varentry{i},'visible','off');
    set(gParams.ui.varname(i),'visible','off');
    set(gParams.ui.incdec{i}{1},'visible','off');
    set(gParams.ui.incdec{i}{2},'visible','off');
  end
end

%--------------------------------- make Help/Ok/Cancel buttons-----------------------------------
gParams.callback = [];
makeOkButton = 1;
makeCancelButton = modal; %only put a cancel if dialog box is modal (or if a custom cancel callback is passed)

if modal==0 && fieldIsNotDefined(gParams,'okCallback') %if the dialog is non-modal, ok just closes it, unless custom ok/cancel callbackz are passed)
  gParams.okCallback = @closeHandler;
  okString = 'Close';
else
  okString = 'OK';
end
% see if this has a callback, in which case we don't
% need to make ok/cancel buttons
if ~isempty(callback)
  gParams.callback = callback;
  % if another argument is specified that should
  % be sent as an argument to the callback function
  if ~isempty(callbackArg)
    gParams.callbackArg = callbackArg;
  end
  params = gParams.fignum;
  params2 = mrParamsGet(vars);
  % if another argument is specified then put up 
  % an ok button with the callback
  if ~isempty(okCallback)
    okString = 'OK';
%   else
%     makeOkButton = 0;
  end
  % if a final argument is specified then put up 
  % a cancel button with the callback
  if ~isempty(cancelCallback)
    makeCancelButton = 1;
  else
    makeCancelButton = 0;
  end
else
  gParams.callback = [];
end

% position of buttons
totalColWidth = 1/dParams.multiCols;
thisButtonWidth = min(100/figWidth,totalColWidth/3);
bottomMargin = uiParams.topMargin/figHeight;
thisButtonHeight = uiParams.buttonHeight/figHeight;
intervalBetweenButtons = (totalColWidth - thisButtonWidth*3)/(4);
leftPosition = (dParams.multiCols - 1)/dParams.multiCols + intervalBetweenButtons;

%Help Button
gParams.fignum(2) = figure('visible','off');
set(gParams.fignum(2),'userdata',0); %this is just to tell helpHandler if the help figure has been drawn or not
gParams.figlocstr{2} = sprintf('mrParamsDialogHelp_%s',fixBadChars(titleStr));
gParams.helpButton = uicontrol(gParams.fignum(1),'Style','pushbutton','Callback',{@helpHandler,gParams.fignum(2),uiParams},'String','Show help',...
  'Position',[leftPosition bottomMargin thisButtonWidth thisButtonHeight],...
  'FontSize',uiParams.fontsize,'FontName',uiParams.fontname);

%Cancel Button
if makeCancelButton
  gParams.cancelButton = uicontrol(gParams.fignum(1),'Style','pushbutton','Callback',@cancelHandler,'String','Cancel',...
  'Position',[leftPosition+(intervalBetweenButtons+thisButtonWidth) bottomMargin thisButtonWidth thisButtonHeight],...
  'FontSize',uiParams.fontsize,'FontName',uiParams.fontname);
end

%Ok Button
if makeOkButton
  gParams.okButton = uicontrol(gParams.fignum(1),'Style','pushbutton','Callback',@okHandler,'String',okString,...
    'Position',[leftPosition+(intervalBetweenButtons+thisButtonWidth)*2 bottomMargin thisButtonWidth thisButtonHeight],...
    'FontSize',uiParams.fontsize,'FontName',uiParams.fontname);
end

%if non-modal, quit here
if ~modal,return,end


% set the input control to the first field that is editable
focusSet = 0;
if isfield(gParams,'ui') && isfield(gParams.ui,'singleEntry')
  % set the first editable field to have the keyboard focus
  for i = 1:length(gParams.varinfo)
    if gParams.varinfo{i}.editable
      % check to see if it is an array of handles
      if isequal(size(gParams.ui.varentry{i}),[1 1])
	H = gParams.ui.varentry{i};
      else
	H = gParams.ui.varentry{i}(1);
      end
      % confirm that we have a handle
      if ishandle(H)
	uicontrol(H);
	focusSet = 1;
	break;
      end
    end
  end
end

% if focus has not been set, then set the focus to the figure
% so the keyboard handler will let you Esc to cancel / enter to ok
if ~focusSet
  figure(gParams.fignum(1));
end

% set keyboard function
if ieNotDefined('ignoreKeys')
  set(gParams.fignum(1),'KeyPressFcn',@mrParamsKeyPressFcn);
end

% wait for user to hit ok or cancel (which sets uiresume)
uiwait;

% this can happen if we have opened up (and closed a second
% mrParamsDialog--this should be removed when this function
% is fixed to be able to run multiple simultaneous mrParamsDialogs;
if ieNotDefined('gParams'),params=[];params2=[];return,end

% check return value
switch(gParams.ok)
  case 1
    params = mrParamsGet(vars);
  case 0     % if cancel has been pressed return empty
    params = [];
  case -1 %if cancel button on top of window has been pressed, return 0*1 empty
    params = zeros(0,1);
end
params2 = [];

closeHandler;


%%%%%%%%%%%%%%%%%%%%
% callback for button handler
%%%%%%%%%%%%%%%%%%%%
function buttonHandler(varnum,entryRow,entryCol,incdec)

global gParams;


%get the value
if strcmp(gParams.varinfo{varnum}.type,'checkbox')
  val = get(gParams.ui.varentry{varnum}(entryRow,entryCol),'Value');
elseif strcmp(gParams.varinfo{varnum}.type,'popupmenu')
  val = [];
  if isfield(gParams.varinfo{varnum},'controls')
    % get the value from the list of values
    val = get(gParams.ui.varentry{varnum}(entryRow,entryCol),'Value');
    val = gParams.varinfo{varnum}.value{val};
    if ischar(val),val=mrStr2num(val);end
  end
else
  % get the value of the text field
  val = get(gParams.ui.varentry{varnum}(entryRow,entryCol),'string');
  if ~any(strcmp(gParams.varinfo{varnum}.type,{'string','stringarray'}))
    % convert to number
    val = mrStr2num(val);
  end
end

isPushButton = strcmp(gParams.varinfo{varnum}.type,'pushbutton');
% if this is a push button and there is no callback
if isPushButton && ~isfield(gParams.varinfo{varnum},'callback')
    disp(sprintf('(mrParamsDialog) Pushbutton %s does not have a callback',gParams.varinfo{varnum}.name));
    return
%if this is a pushbutton or the value is returned by a callback
elseif isfield(gParams.varinfo{varnum},'callback') && ...
        isPushButton || gParams.varinfo{varnum}.passCallbackOutput
  args = {};%getVars = 0;
  % if it wants optional arguments, pass that
  if isfield(gParams.varinfo{varnum},'callbackArg')
    args{end+1} = gParams.varinfo{varnum}.callbackArg;
  end
  % if the function wants the current parameter settings, pass that
  if isfield(gParams.varinfo{varnum},'passParams') && (gParams.varinfo{varnum}.passParams == 1)
    args{end+1} = mrParamsGet(gParams.vars);
  end
  % if the function wants the entry value
  if isfield(gParams.varinfo{varnum},'passValue') && (gParams.varinfo{varnum}.passValue == 1)
    args{end+1} = val;
  end

  %call the function
  if gParams.varinfo{varnum}.passCallbackOutput
    val = callbackEval(gParams.varinfo{varnum}.callback,args);
    if isPushButton
      gParams.varinfo{varnum}.value = val;
      return
    end
  else
    callbackEval(gParams.varinfo{varnum}.callback,args);
    return; %if no value is returned, we're done
  end
end

% if this is supposed to be a number, then make sure it is.
if ~any(strcmp(gParams.varinfo{varnum}.type,{'string','stringarray'}))
  % check for incdec (this is for buttons that increment or decrement
  % the values, if one of these was passed, we will have by how much
  % we need to increment or decrement the value
  if isfield(gParams.varinfo{varnum},'incdec') && exist('incdec','var')
    val = val+incdec;
  end
  % check if number needs to be round
  if isfield(gParams.varinfo{varnum},'round') && gParams.varinfo{varnum}.round
    roundVal = gParams.varinfo{varnum}.round;
    if (round(roundVal) == roundVal) && (roundVal > 1)
      % round to a significant number of digits
      val = round(val*(10^(roundVal-1)))/(10^(roundVal-1));
    else
      % round if set to 1
      val = round(val);
    end
  end
  % check for minmax violations
  if isfield(gParams.varinfo{varnum},'minmax')
    % check minimum value (second check is to make sure it has exceeded minimum by more
    % than a tiny round-off problem)
    if (val < gParams.varinfo{varnum}.minmax(1)) && ((gParams.varinfo{varnum}.minmax(1)-val) > 10e-12)
      disp(sprintf('(mrParamsDialog) Value %f lower than minimum %f',val,gParams.varinfo{varnum}.minmax(1)));
      val = [];
    % check maximum value (second check is to make sure it has exceeded minimum by more
    % than a tiny round-off problem)
    elseif (val > gParams.varinfo{varnum}.minmax(2)) && ((val - gParams.varinfo{varnum}.minmax(1)) > 10e-12)
      disp(sprintf('(mrParamsDialog) Value %f greater than maximum %f',val,gParams.varinfo{varnum}.minmax(2)));
      val = [];
    end
  end
  % if the variable is empty, then set it back to default, this happens
  % if we got an invalid number entry
  if isempty(val)
    if ~strcmp(gParams.varinfo{varnum}.type,'popupmenu')
      if ~isempty(gParams.varinfo{varnum}.value)
        set(gParams.ui.varentry{varnum}(entryRow,entryCol),'string',gParams.varinfo{varnum}.value(entryRow,entryCol));
      end
      return  %no need to go to callback if the value hasn't changed, so return
    end
    % otherwise remember this string as the default
  else
    if ~any(strcmp(gParams.varinfo{varnum}.type,{'popupmenu','checkbox'}))
%       if ismember(gParams.varinfo{varnum}.type,{'array','numeric'})
        gParams.varinfo{varnum}.value(entryRow,entryCol)=val;
%       else
%         gParams.varinfo{varnum}.value(entryRow,entryCol) = thisNum2str(val);
%       end
      set(gParams.ui.varentry{varnum}(entryRow,entryCol),'string',thisNum2str(val));
    end
    % now check to see if this variable controls another one
    if isfield(gParams.varinfo{varnum},'controls')
      % go through all the controlled values
      for i = gParams.varinfo{varnum}.controls
        % enable or disable all the controlled fields
        if (val == 0)
	  % if set to not, then when this is off the controls is on
	  if gParams.varinfo{varnum}.controlsNot
	    set(gParams.ui.varentry{i},'Enable','on');
	  else
	    set(gParams.ui.varentry{i},'Enable','off');
	  end
        else
	  % if set to not, then when this is off the controls is on
	  if gParams.varinfo{varnum}.controlsNot
	    set(gParams.ui.varentry{i},'Enable','off');
	  else
	    set(gParams.ui.varentry{i},'Enable','on');
	  end
        end
        % now store the value currently being displayed
        if gParams.varinfo{i}.oldControlVal
          % get the current value according to the GUI
          currentValue = get(gParams.ui.varentry{i},'String');
	  % for some types you need to parse things properly to
	  % get the current value
          if strcmp(gParams.varinfo{i}.type,'popupmenu')
	    valueNum = get(gParams.ui.varentry{i},'Value');
	    % valueNum should not be 0, but sometimes it can get in that
	    % way if you click funny on the popupmenu and it doesn't select
	    % anything, this is a fix from getting lots of debug messages
	    % when that happens
	    if valueNum > 0
	      currentValue = putOnTopOfList(currentValue{valueNum},currentValue);
	    else
	      currentValue = putOnTopOfList(currentValue{1},currentValue);
	    end
          elseif strcmp(gParams.varinfo{i}.type,'checkbox')
            currentValue = get(gParams.ui.varentry{i},'Value');
          elseif strcmp(gParams.varinfo{i}.type,'pushbutton')
	    currentValue = gParams.varinfo{i}.value;
          end
          % and save it
	  if strcmp(gParams.varinfo{i}.type,'array')
	    for k = 1:length(currentValue)
	      currentNumericValue(k) = str2num(currentValue{k});
	    end
	    currentNumericValue = reshape(currentNumericValue,size(currentValue));
	    gParams.varinfo{i}.allValues{gParams.varinfo{i}.oldControlVal}=currentNumericValue;
	  else
	    gParams.varinfo{i}.allValues{gParams.varinfo{i}.oldControlVal}=currentValue;
	  end
        end
        % switch to new value
        if (val >=1) && (val <= length(gParams.varinfo{i}.allValues))
          gParams.varinfo{i}.value = gParams.varinfo{i}.allValues{val};
	  % if this is an array, we have to set each individual array item
	  if strcmp(gParams.varinfo{i}.type,'array')
	    for k = 1:length(gParams.varinfo{i}.allValues{val})
	      set(gParams.ui.varentry{i}(k),'String',gParams.varinfo{i}.allValues{val}(k));
	    end
	  else
	    if ~any(strcmp(gParams.varinfo{i}.type,{'checkbox','pushbutton'})) 
	      set(gParams.ui.varentry{i},'String',gParams.varinfo{i}.allValues{val});
	    end
	  end
	  % some more things to set for these types
	  if strcmp(gParams.varinfo{i}.type,'popupmenu')
	    set(gParams.ui.varentry{i},'Value',1);
	  elseif strcmp(gParams.varinfo{i}.type,'checkbox')
	    if ischar(gParams.varinfo{i}.value)
	      set(gParams.ui.varentry{i},'Value',str2num(gParams.varinfo{i}.value));
	    else
	      set(gParams.ui.varentry{i},'Value',gParams.varinfo{i}.value);
	    end
    end
          gParams.varinfo{i}.oldControlVal = val;
        end
      end
    end
  end
  % if the field has incdec, see how they should be grayed or not
  if isfield(gParams.varinfo{varnum},'incdec') && isfield(gParams.varinfo{varnum},'minmax')
    enableArrows(val,varnum,entryRow,entryCol)
  end
end
% update params
if isfield(gParams, 'callback')
  if ~isempty(gParams.callback)
    gParams.params = mrParamsGet(gParams.vars);
    if isfield(gParams,'callbackArg')
      callbackEval(gParams.callback,{gParams.params},{gParams.callbackArg}); 
    else
      callbackEval(gParams.callback,{gParams.params});
    end
  end
end

% handle callbacks for non-push buttons
if ~ieNotDefined('gParams') && isfield(gParams.varinfo{varnum},'callback')...
    && ~gParams.varinfo{varnum}.passCallbackOutput && ~strcmp(gParams.varinfo{varnum}.type,'pushbutton')
  callbackArgs={};
  if isfield(gParams.varinfo{varnum},'callbackArg')
    callbackArgs{end+1} = gParams.varinfo{varnum}.callbackArg;
  end
  callbackArgs{end+1} = mrParamsGet(gParams.vars);
  callbackEval(gParams.varinfo{varnum}.callback,callbackArgs);
  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% turn on or off incdec arrows depending on minmax
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function enableArrows(val,varnum,entryRow,entryCol)

global gParams;

if isfield(gParams.varinfo{varnum},'incdec') && isfield(gParams.varinfo{varnum},'minmax')
  if isnumeric(val)
    % turn on or off dec arrow
    if (val+gParams.varinfo{varnum}.incdec(1)) < gParams.varinfo{varnum}.minmax(1)
      set(gParams.ui.incdec{varnum}{1}(entryRow,entryCol),'Enable','off');
    else
      set(gParams.ui.incdec{varnum}{1}(entryRow,entryCol),'Enable','on');
    end
    % turn on or off inc arrow
    if (val+gParams.varinfo{varnum}.incdec(2)) > gParams.varinfo{varnum}.minmax(2)
      set(gParams.ui.incdec{varnum}{2}(entryRow,entryCol),'Enable','off');
    else
      set(gParams.ui.incdec{varnum}{2}(entryRow,entryCol),'Enable','on');
    end
  end
end

%%%%%%%%%%%%%%%%%%%%
% callback for help
%%%%%%%%%%%%%%%%%%%%
function helpHandler(handle,event,fignum,uiParams)

if strcmp(get(fignum,'visible'),'on')
  set(fignum,'visible','off');
  set(handle,'string','Show Help');
else
  set(fignum,'visible','on');
  set(handle,'string','Hide Help');
end

if get(fignum,'userdata')

else  
  global gParams

  % display to text buffer
  for i = 1:length(gParams.varinfo)
    disp(sprintf('%s: %s\n',gParams.varinfo{i}.name,gParams.varinfo{i}.description));
  end
  
  % turn off menu/title etc.
  set(fignum,'MenuBar','none');
  set(fignum,'NumberTitle','off');
  set(fignum,'Name','Parameter help');
  set(fignum,'closeRequestFcn',@helpcloseHandler);

  set(fignum,'defaultUicontrolUnits','pixels'); 
  
  dParams.entryWidth = zeros(1,length(gParams.varinfo));
  dParams.entryValue = zeros(1,length(gParams.varinfo));
  dParams.entryString = repmat({{''}},1,length(gParams.varinfo));
  dParams.entryStyle =  repmat({'text'},1,length(gParams.varinfo));
  dParams.incdec = zeros(length(gParams.vars),2);
  dParams.testString = repmat({''},1,length(gParams.varinfo));
  dParams.entryNumCols = ones(1,length(gParams.varinfo));
  dParams.entryNumRows = ones(1,length(gParams.varinfo));
  dParams.numLines = ones(1,length(gParams.varinfo));
  for i = 1:length(gParams.varinfo)
    if ~gParams.varinfo{i}.visible %no need to display the help is parameter is not visible
      dParams.numLines(i)=0;
    end
    dParams.entryString{i} = {[' ' gParams.varinfo{i}.description]}; %add 1 space on the left
    dParams.testString(i) = dParams.entryString{i};
  end

  %gParams.fignum(2) = fignum;
  
  %compute figure dimensions based on number of rows and colums
  [figurePosition,dParams, uiParams] = optimizeFigure(fignum,gParams.figlocstr{2},dParams,uiParams);
  
  %set the all the entry widths to the max 
  dParams.entryWidth(:)=dParams.allEntriesWidth;
  
  %set control dimensions to normalized so that the figure resizes
  set(fignum,'defaultUicontrolUnits','normalized'); 
  % put up the info
  for i = 1:length(gParams.varinfo)
    if gParams.varinfo{i}.visible %no need to display the help is parameter is not visible
      makeUIcontrol(i,fignum,dParams,uiParams,'varname');
      set(makeUIcontrol(i,fignum,dParams,uiParams,'varentry'),'HorizontalAlignment','Left');
    end
  end

  % make close button
  uicontrol(fignum,'Style','pushbutton','Callback',@helpcloseHandler,'String','Close',...
    'Position',getUIControlPos(fignum,dParams,uiParams,dParams.figrows,1,dParams.multiCols,2,dParams.allEntriesWidth,1),...
    'FontSize',uiParams.fontsize,'FontName',uiParams.fontname);
  
  set(fignum,'userdata',1)
end

%%%%%%%%%%%%%%%%%%%%
% callback for helpclose
%%%%%%%%%%%%%%%%%%%%
function helpcloseHandler(varargin)

global gParams;

set(gParams.fignum(2),'visible','off')
if ishandle(gParams.helpButton)
  set(gParams.helpButton,'string','Show Help');
end

%%%%%%%%%%%%%%%%%%%%
% callback for close
%%%%%%%%%%%%%%%%%%%%
function closeHandler(varargin)

global gParams;
if isempty(gParams),return,end

if isfield(gParams,'fignum') 
  if isfield(gParams,'figlocstr')
  % save figure locations .mrDefaults
    for iFig = 1:length(gParams.fignum)
      if ishandle(gParams.fignum(iFig)) && (length(gParams.figlocstr)>=iFig)
	mrSetFigLoc(fixBadChars(gParams.figlocstr{iFig}),get(gParams.fignum(iFig),'Position'));
      end
    end
  end
  % close figure if it exists
  if  ishghandle(gParams.fignum)
    delete(gParams.fignum);
  end
end
saveMrDefaults;

clear global gParams;
drawnow


%%%%%%%%%%%%%%%%%%%%
% callback for ok
%%%%%%%%%%%%%%%%%%%%
function okHandler(varargin)

global gParams;
gParams.ok = 1;
if isfield(gParams,'okCallback')
  if isfield(gParams,'callbackArg')
    callbackEval(gParams.okCallback,[],gParams.callbackArg);
  else
    callbackEval(gParams.okCallback);
  end
  closeHandler;
else
  uiresume;
end

%%%%%%%%%%%%%%%%%%%%
% callback for cancel
%%%%%%%%%%%%%%%%%%%%
function cancelHandler(varargin)

global gParams;
gParams.ok = 0;
if isfield(gParams,'cancelCallback')
  if isfield(gParams,'callbackArg')
    callbackEval(gParams.cancelCallback,[],gParams.callbackArg);
  else
    callbackEval(gParams.cancelCallback);
  end
  closeHandler;
else
  uiresume;
end

%%%%%%%%%%%%%%%%%%%%
% callback for close button on top of window
%%%%%%%%%%%%%%%%%%%%
function altCloseHandler(varargin)

global gParams;
gParams.ok = -1;
uiresume;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% makeUIcontrol makes an uicontrol of any type %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [hEntry,hMinus,hPlus] = makeUIcontrol(varnum,fignum,dParams,uiParams,columnType)

hMinus = [];
hPlus = [];

multiCol = find(dParams.startMultiCol<=varnum,1,'last');
numRows = dParams.numLines.*dParams.entryNumRows;
switch(columnType)
  case 'varname'
    colnum=1;
    entryWidth = uiParams.varNameWidth;
    fieldnums = 1;
    rownums = sum(numRows(dParams.startMultiCol(multiCol):varnum-1))+1;
    numLines = 	dParams.numLines(varnum)*dParams.entryNumRows(varnum);
    entryString = uiParams.varName(varnum);
    style = 'text'; 
    hAlignment = 'right';
    incdec = [0 0];
    
  case 'varentry'
    colnum =2;
    entryWidth = dParams.entryWidth(varnum);
    fieldnums =1:dParams.entryNumCols(varnum);
    rownums = sum(numRows(dParams.startMultiCol(multiCol):varnum-1))+(1:dParams.entryNumRows(varnum));
    numLines = 	dParams.numLines(varnum);
    entryString = dParams.entryString{varnum};
    style = dParams.entryStyle{varnum};
    hAlignment = 'center';
    incdec = dParams.incdec(varnum,:);
    
end
if ~numLines %if numLines==0, that means the control is invisble
  numLines =1;  %set it to 1 to avoid an error
end


for i=1:length(rownums)
  for j=fieldnums
    
    uiPosition = getUIControlPos(fignum,dParams,uiParams,rownums(i),numLines,multiCol,colnum,entryWidth,j);
    
    if isnumeric(entryString{i,j})
      %convert numerical values into string using different precision if they're integer or decimal
      entryString{i,j} = thisNum2str(entryString{i,j}); 
    end
    hEntry(i,j) = uicontrol(fignum,...
    'Style',style,...
    'Callback',sprintf('mrParamsDialog(%f,%f,%f)',varnum,i,j),...  %callback has no effect if textbox
    'String',entryString{i,j},...
    'Value',dParams.entryValue(varnum),...
    'Position',uiPosition,...
    'HorizontalAlignment',hAlignment,...
    'FontSize',uiParams.fontsize,'FontName',uiParams.fontname);
  
    if any(incdec)
      % make callback string
      deccallback = sprintf('mrParamsDialog(%f,%f,%f,%f)',varnum,i,j,incdec(1));
      inccallback = sprintf('mrParamsDialog(%f,%f,%f,%f)',varnum,i,j,incdec(2));
      
      figurePosition = get(fignum,'position');
      incdecMargin = uiParams.incdecMargin/figurePosition(3);
      entryPosition =get(hEntry(i,j),'position');
      decPosition = entryPosition;
      incPosition = entryPosition;
      %compute new positions
          
      switch(dParams.incdecType{varnum})
        case 'plusMinus'
          incdecWidth = min(uiParams.maxIncdecButtonWidth+incdecMargin,entryWidth/2)/figurePosition(3)-incdecMargin;
          entryPosition(3) = entryPosition(3)-(incdecWidth+incdecMargin);
          decPosition(1) = entryPosition(1)+entryPosition(3)+incdecMargin;
          if ismac 
            incPosition(2) = incPosition(2)+incPosition(4)*.5;
            decPosition(4) = decPosition(4)*.45;
            incPosition(4) = incPosition(4)*.45;
          else %on Windows, the minus sign is too low, so I'm making the button taller so that we can see it
            %also the button width is smaller
            incPosition(2) = incPosition(2)+incPosition(4)*.45;
            decPosition(4) = decPosition(4)*.75;
            incPosition(4) = incPosition(4)*.55;
          end
          incString = '+';
          decString = '-';

        case 'arrows'
          incdecWidth = min(uiParams.maxIncdecButtonWidth+incdecMargin,entryWidth/3)/figurePosition(3)-incdecMargin;
          entryPosition(1) = entryPosition(1)+incdecWidth+incdecMargin;
          entryPosition(3) = entryPosition(3)-2*(incdecWidth+incdecMargin);
          incString = '>';
          decString = '<';

      end
      incPosition(1) = entryPosition(1)+entryPosition(3)+incdecMargin;
      decPosition(3) = incdecWidth;
      incPosition(3) = incdecWidth;

      % if width goes below this amount, then skip the +/- buttons
      minUISize = 0.001;
      if (incdecWidth < minUISize)
	% Resize entry to not include incdec
	entryPosition(3) = entryPosition(3)+(incdecWidth+incdecMargin);
	set(hEntry(i,j),'position',entryPosition);
	% remove the incdec field from gParams
	global gParams;
	if isfield(gParams.varinfo{varnum},'incdec')
	  gParams.varinfo{varnum} = rmfield(gParams.varinfo{varnum},'incdec');
	end
      else
	%correct position of entry field
	set(hEntry(i,j),'position',entryPosition);
	%make decrement and increment buttons
	hMinus(i,j) = uicontrol(fignum,'Style','pushbutton','Callback',deccallback,'String',decString,...
				'Position',decPosition,'FontSize',uiParams.fontsize,'FontName',uiParams.fontname);
	hPlus(i,j) = uicontrol(fignum,'Style','pushbutton','Callback',inccallback,'String',incString,...
			       'Position',incPosition,'FontSize',uiParams.fontsize,'FontName',uiParams.fontname);
      end
    end
  end
end

if colnum==2 && strcmp(dParams.entryStyle{varnum},'checkbox') %on windows, the backgroud of checkboxes is colored
  set(hEntry,'BackgroundColor',get(gcf,'color')); %even without string, which is ugly
end
if strcmp(dParams.entryStyle{varnum},'text') && isempty(entryString{1}) %if it's an empty textbox
  set(hEntry,'visible','off');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% getUIControlPos returns a location for a uicontrol %
%   dealing with multicolumns and margins            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pos = getUIControlPos(fignum,dParams,uiParams,rownum,numLines,multiCol,colnum,entryWidth,fieldNum)

numcols = 2;
% always make sure that the last row end up on the
% last row even if we have multiple columns
if multiCol == dParams.multiCols && (rownum+numLines-1) == dParams.figrows
    rownum = dParams.figrows-numLines+1;
end
colnum = numcols*(multiCol-1)+colnum;

% get figure position
figurePosition = get(fignum,'Position');

% set the horizontal position and width for the button
if colnum - (multiCol-1)*numcols == 1
  pos(1) = uiParams.leftMargin + ... %position
            (multiCol - 1) * (uiParams.varNameWidth + uiParams.margin) + ...
            (multiCol - 1) * (dParams.allEntriesWidth+uiParams.margin) + ...
            uiParams.margin; 
else
  pos(1) = uiParams.leftMargin + ...
    multiCol*(uiParams.margin + uiParams.varNameWidth) + ...
    (colnum-multiCol-1) * (dParams.allEntriesWidth+uiParams.margin) + ...
    (fieldNum-1)*(entryWidth+uiParams.margin)+...
    uiParams.margin;
end
pos(3) = entryWidth; 

% set the vertical position and height for the button
pos(4) = uiParams.buttonHeight*numLines+uiParams.margin*(numLines-1);
pos(2) = figurePosition(4)-pos(4)-uiParams.topMargin - (uiParams.buttonHeight+uiParams.margin)*(rownum-1);

%normalize position
pos([1 3]) = pos([1 3])/figurePosition(3);
pos([2 4]) = pos([2 4])/figurePosition(4);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% optimizeFigure optimizes the number of rows and columns as well as the dimensions of the figure %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [figurePosition,dParams,uiParams] = optimizeFigure(fignum,figLocStr,dParams,uiParams)

%compute figure dimensions based on number of rows and colums
figurePosition = mrGetFigLoc(fixBadChars(figLocStr));
if isempty(figurePosition)
  figurePosition = get(fignum,'Position');
end

maxEntryNumCols = max(dParams.entryNumCols);

%-------------------optimize figure dimensions 
monitorPositions = getMonitorPositions;
[whichMonitor,figurePosition]=getMonitorNumber(figurePosition,monitorPositions);
screenSize = monitorPositions(whichMonitor,:); % find which monitor the figure is displayed in

dParams.multiCols=0;                             %
actualMultiCols = 0;                             %these values are set only to pass the first test
figHeight = uiParams.maxFigHeightWidthRatio+1;   %
figWidth = 1;                                    %
triedReduceEntryWidth = 0;
minFontSize = 6;
%while one of the dimensions is larger than the screen or the height/width is over the threshold, resize
while figHeight/figWidth>uiParams.maxFigHeightWidthRatio || figWidth>screenSize(1,3) || figHeight>screenSize(1,4)

  % test is figure respects constraints, if not, change something
  if figWidth<screenSize(1,3) && dParams.multiCols<length(dParams.numLines)+1
  % first try adding columns (no more than the number of entries), but only if the screen width has not been reached
    dParams.multiCols = dParams.multiCols+1;
  elseif figWidth>screenSize(1,3) && ~triedReduceEntryWidth
  %then try reducing  the max entry width (only if figure width is large than screen width)
    %compute the max entries width as what's left to the screen width when you remove varnames and margins and divide by multicols
    uiParams.maxEntriesWidth = floor((screenSize(1,3)- 2*uiParams.leftMargin ...
                              - actualMultiCols*(uiParams.varNameWidth+uiParams.margin)...
                              - (actualMultiCols-1)*uiParams.margin)...
                              / actualMultiCols);
    uiParams.maxSingleFieldWidth = uiParams.maxEntriesWidth / maxEntryNumCols;
    %set the min entry width to 0
    uiParams.minEntriesWidth = 0;
    triedReduceEntryWidth=1;
    
  elseif uiParams.fontsize>minFontSize
  %finally try reduce the fontsize
    uiParams.fontsize = uiParams.fontsize-1;
  else
    %if everything has been tried, break out of the while loop
    break;
  end

  %-------------------compute the uicontrol and figure dimensions
  %%%%%%%%%%%% get varname string width for fields that might wrap
  uiParams.varNameWidth =0;
  for i = 1:length(uiParams.varName)  
    h = uicontrol(fignum,'Style','text','String',uiParams.varName{i},'FontSize',uiParams.fontsize,'FontName',uiParams.fontname);
    thisExtent = get(h,'extent');
    uiParams.varNameWidth = max(uiParams.minVarNameWidth,max(thisExtent(3),uiParams.varNameWidth));
    delete(h);
  end
  
  %%%%%%%%%%%% get string width for fields that might wrap
  for i = 1:length(dParams.testString)  
    if ~isempty(dParams.testString{i}) && dParams.numLines(i)~=0
      %compute number of lines using string width if it's gonna be displayed using a text box, a popupmenu or a pushbutton
      h = uicontrol(fignum,'Style',dParams.entryStyle{i},'String',dParams.testString{i},'FontSize',uiParams.fontsize,'FontName',uiParams.fontname);
      thisExtent = get(h,'extent');
      dParams.entryWidth(i) = thisExtent(3)+20; %we need to allow some space for the button features
      delete(h);
    end
  end
  
  %%%%%%%%%%%% get entry height
  if ieNotDefined('thisExtent')
    h = uicontrol(fignum,'Style','Text','String','X','FontSize',uiParams.fontsize,'FontName',uiParams.fontname);
    thisExtent = get(h,'extent');
    delete(h);
  end
  uiParams.buttonHeight = thisExtent(4);
  %For edit boxes and buttons on MACs, this height will be too small because of their large borders
  if ismac 
    uiParams.buttonHeight = uiParams.buttonHeight*1.25;
  end
  % global mrDEFAULTS;                % The height of the button used to be dependent on the version of matlab             
  % mver = matlabVersionNumber;       % in addition ot the computer type. not sure this is useful anymore
  % if strcmp(computer,'MACI') || strcmp(computer,'MACI64') || (mver > 7.4)
  %   ...
  
  %%%%%%%%%%%% compute the total entry width
  %the total field width is whatever field has the largest width, within the min and max parameters
  dParams.allEntriesWidth = max(max(dParams.entryWidth),min(maxEntryNumCols*uiParams.maxSingleFieldWidth,uiParams.maxEntriesWidth));
  dParams.allEntriesWidth = max(uiParams.minEntriesWidth,min(uiParams.maxEntriesWidth,dParams.allEntriesWidth));

  %%%%%%%%%%%% get  number of lines for fields that might wrap
  for i = 1:length(dParams.entryStyle)  
    if ~isempty(dParams.testString{i}) && dParams.numLines(i)~=0 && ~strcmp(dParams.entryStyle{i},'popupmenu')
      dParams.numLines(i) = ceil(ceil(dParams.entryWidth(i)/dParams.allEntriesWidth)*uiParams.lineHeightRatio);
    end
  end

  %%%%%%%%%%%% compute total number of entries per multicols
  numRows = [dParams.numLines.*dParams.entryNumRows 1]; %we add one for the help/ok/cancel buttons
  %compute new number of rows per columns, but make sure we're not cutting an entry
  dParams.figrows = max(max(numRows),ceil(sum(numRows)/dParams.multiCols)) - 1; %we remove one just because of the order of things in the while loop
  %now check if any multirow entry is split between two columns and if we have enough rows per column
  keepAddingRows=1; %this is just to enter the while loop
  while keepAddingRows   %while an entry is split between two columns or 
    %there are more rows than number of columns times number of rows per column
    dParams.figrows = dParams.figrows+1;  % we try adding one row per column
    numRowsLeft = numRows;
    dParams.startMultiCol = zeros(1,dParams.multiCols);
    endMultiCol = 0;
    cutsEntry=0;
    iCol=0;
    while iCol<dParams.multiCols && ~isempty(numRowsLeft) %and we check for each column if an entry is split. 
      iCol=iCol+1;  % we do this in a while rather than a for loop because there are cases where adding a row removes more than one column at once
                    % in which case the last column can be empty
      %this is the entry index of the start of the current column (the first entry 
      % whose cumulated number of rows is less than the number of rows per column)
      dParams.startMultiCol(iCol) = endMultiCol+find(cumsum(numRowsLeft)<=dParams.figrows,1,'first');
      % similarly, this is the last entry whose cumulated number of rows is less than the number of rows per column
      endMultiCol = endMultiCol+find(cumsum(numRowsLeft)<=dParams.figrows,1,'last');
      %these are the number of rows of all the entries attributed to the current columns
      thisNumRows = numRows(dParams.startMultiCol(iCol):endMultiCol);
      % and the rows left for next columns
      numRowsLeft = numRows(endMultiCol+1:end);
      %an entry is split between two columns if it was split for any previous column
      %or the cumulated number of rows per entry (from the start of the current column)
      %does not include the total number of rows per column and the column is full 
      %(if the total number of rows is less than the desired number, then an entry cannot be split over the next column)
      cutsEntry = cutsEntry || ~(ismember(dParams.figrows,cumsum(thisNumRows)) || dParams.figrows>=sum(thisNumRows));
    end
    % We remember the actual number of column needed, in case we got out of the while loop before reaching the last column
    actualMultiCols=iCol; 
    %check that there are not entries left (if there are, we need to increase the number of rows even if there are no cut entries)
    keepAddingRows = cutsEntry || ~isempty(numRowsLeft);
  end

  %%%%%%%%%%%% compute figure dimensions
  figHeight = 2*uiParams.topMargin+dParams.figrows*uiParams.buttonHeight+(dParams.figrows-1)*uiParams.margin;
  figWidth = 2*uiParams.leftMargin...
              + (actualMultiCols- 1)*uiParams.margin...
              + actualMultiCols*(uiParams.varNameWidth+dParams.allEntriesWidth+uiParams.margin);

end

dParams.multiCols = actualMultiCols;

% set the figure position
figurePosition(4) = figHeight;
figurePosition(3) = figWidth;
%make sure the figure is not outside the screen
figurePosition(1) = min(figurePosition(1),sum(screenSize([1 3]))-1-figWidth);
figurePosition(2) = min(figurePosition(2),sum(screenSize([2 4]))-1-figHeight);

set(fignum,'Position',figurePosition);

%replace non-set widths by the max width
dParams.entryWidth(dParams.entryWidth<0)= dParams.allEntriesWidth;

%modified num2str to increase the number of decimals for reals
function value = thisNum2str(value)

  if rem(value,1)~=0
    value = num2str(value,'%.6f');
  else
    value = num2str(value);
  end





