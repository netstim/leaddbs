% mlrAdjustGUI.m
%
%      usage: mlrAdjustGUI(v,command,varargin)
%         by: justin gardner
%       date: 10/31/10
%    purpose: Adjusts the MLR GUI to allow plug-in code
%             to change the behavior of the GUI. Note
%             that this serves a different function then
%             mlrGuiSet which is used to set user interface
%             items during normal operations (like graying
%             out a menu item, or changing the value of a control).
%             This is intended to be used only when MLR is loaded
%             to adjust the GUI to suit the needs of different sites.
%
%             Every UI interface item on the MLR window
%             can be specified by either a tag (a string that
%             identifies it), or by the label for a menu item.
%             e.g., the menu item under ROI called Show
%             is /ROI/Show. or e.g., the cortical depth slider
%             has the tag corticalDepth. These identifiers are
%             what is called the "controlName" below and you
%             can get a full list of all possible controlNames 
%             by doing:
%             controlList = mlrAdjustGUI(getMLRView,'list');
%
%             To set a property of a control: 
%             mlrAdjustGUI(v,'set',controlName or tag,propertyName,propertyValue);
%      e.g.:  mlrAdjustGUI(getMLRView,'set','baseGammaSlider','Visible','off');
% 
%             Similarly, to set the callback for a control
%      e.g.:  mlrAdjustGUI(getMLRView,'set','baseGammaSlider','Callback',@testCallback);
%             Where testCallback is a function that takes two arguments and
%             usually begins with a snippet of code that gets the view:
%             function testCallback(hObject,eventdata)
%             v = viewGet(getfield(guidata(hObject),'viewNum'),'view');
%
%             To add a new menu item:
%             mlrAdjustGUI(v,'add','menu',menuName,menuLocation,propertyName1,propertyValue1,...);
%      e.g.:  mlrAdjustGUI(getMLRView,'add','menu','Plugin','Plots','Callback',@testCallback,'Separator','on');
%             The above will add the menu item Plugin after the menu identified
%             as Plots. If you wanted instead to put it at the top
%             of the Plots menu, then set the menuLocation to /Plots/
%             An alternative way is to use a tag to add a menu, and specify the name and tag property 
%             mlrAdjustGUI(v,'add','menu',menuTag,menuLocation,'label',menuName,'tag',menuTag...);
%      e.g.:  mlrAdjustGUI(getMLRView,'add','menu','pluginMenuItem','Plots','label','Plugin','tag','pluginMenuItem');
%             It is easier to make tags unique
% 
%             To move a menu item around
%             mlrAdjustGUI(v,'set',menuName or tag,'location',menuLocation);
%
%             to remove a menu item
%             mlrAdjustGUI(v,'remove','menu',menuName or tag)
%
%             To add an interrogator function as a default one
%             (that shows up in the GUI)
%             mlrAdjustGUI(v,'add','interrogator',interrogatorName)
%      e.g.:  mlrAdjustGUI(getMLRView,'add','interrogator','eventRelatedPlot');
%
%             To add colormap functions which will show up in
%             /Edit/Overlay
%             mlrAdjustGUI(v,'add','colormap',colormapName)
%       e.g.: mlrAdjustGUI(getMLRView,'add','colormap','gray');
%
function retval = mlrAdjustGUI(v,command,varargin)

verbose=false;

% default return empty
retval = [];

% check arguments
if nargin < 2
  help mlrAdjustGUI
  return
end

% get figure
if isempty(v),disp(sprintf('(mlrAdjustGUI) Empty view. Is MLR closed?'));return,end
f = viewGet(v,'fignum');
if isempty(f),disp(sprintf('(mlrAdjustGUI) Passed in view does not have a figure associated with it'));return;end

% get controls and menus
plotAxes = getAxes(f);
uiControls = getUiControls(f);
menuControls = getMenuControls(f);

% do the commanded action
switch command
 case {'set'}
  setItemProperty(varargin,uiControls,menuControls,plotAxes,verbose)
 case {'list'}
  retval = listControlNames(uiControls,menuControls,plotAxes);
 case {'add'}
  switch varargin{1}
    case {'menu'}
     placeMenu({varargin{2:end}},menuControls,'add',verbose,f);
    case {'interrogator','interrogators'}
     addInterrogator(v,varargin{2},verbose);
    case {'colormap','colormaps'}
     addColormap(v,varargin{2},verbose);
    case {'control'}
     retval = addControl(f,{varargin{2:end}},uiControls,v,verbose);
    case {'axes'}
     addAxes(f,{varargin{2:end}},plotAxes,verbose);
    case {'panel'}
     addPanel(f,{varargin{2:end}},v,verbose);
    otherwise
      mrWarnDlg(['(mlrAdjustGUI) Unknow object type ' varargin{1}]);
  end
 case {'remove'}
  switch varargin{1}
    case {'menu'}
     removeMenu(cellArray(varargin{2}),menuControls,verbose);
  end
 case {'get'}
  % return handle, check menu items and ui controls
  retval = getHandle(varargin{1},menuControls,uiControls,plotAxes);
 otherwise
  disp(sprintf('(mlrAdjustGUI) Unknown command: %s',command));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  listControlNames   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%
function retval = listControlNames(uiControls,menuControls,plotAxes)

% list Axes tags
disp('==============    Axes     ============== ');
for i = 1:length(plotAxes)
  disp(sprintf('(%i) tag: %s',i,plotAxes(i).tag));
end

% list uiControl tags
disp('============== UI Controls ============== ');
for i = 1:length(uiControls)
  disp(sprintf('(%i) tag: %s Style: %s',i,uiControls(i).tag,get(uiControls(i).h,'Style')));
end

% list menu labels and tags
disp('============== Menu Items ============== ');
for i = length(menuControls):-1:1
  disp(sprintf('(%i) Menu label: ''%s'' tag: %s',length(menuControls)-i+1,menuControls(i).fullLabel,menuControls(i).tag));
end

retval.uiControls=uiControls;
retval.menuControls=menuControls;
retval.plotAxes=plotAxes;

%%%%%%%%%%%%%%%%%%%
%%%   getHandle  %%
%%%%%%%%%%%%%%%%%%%
function h = getHandle(itemName,controls,controls2,controls3)

h = [];

% check if there is a match in the following fields
searchFields = {'tag','label','fullLabel'};

for i = 1:length(searchFields)
  % check for existence of field
  if isfield(controls,searchFields{i})
    % and if there is a match
    [tf itemNum] = find(strcmp(itemName,{controls.(searchFields{i})}));
    if tf
      h = controls(itemNum).h;
      return
    end
  end
end
  
% if we didn't find anything and we were passed two
% sets of controls, then check second set
if isempty(h) && (nargin >=3)
  h = getHandle(itemName,controls2);
end
%same for the 
if isempty(h) && (nargin ==4)
  h = getHandle(itemName,controls3);
end

% if we didn't find anything and the thing ends with '/' then
% try without the '/'
if isempty(h) && (length(itemName)>0) && (itemName(end) == '/')
  if nargin == 3
    h = getHandle(itemName(1:end-1),controls,controls2);
  else
    h = getHandle(itemName(1:end-1),controls);
  end
end


%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  addInterrogator   %%
%%%%%%%%%%%%%%%%%%%%%%%%%
function addInterrogator(v,interrogatorList,verbose)

% make into a cell array and set in view
interrogatorList = cellArray(interrogatorList);
viewSet(v,'defaultInterrogators',interrogatorList);

% display what we are doing
if verbose, disp(sprintf('(mlrAdjustGUI) Adding default interrogators: %s',cellToCommaDelimited(interrogatorList))); end

%%%%%%%%%%%%%%%%%%%%%
%%%  addColormap   %%
%%%%%%%%%%%%%%%%%%%%%
function addColormap(v,colormapList,verbose)

% make into a cell array and set in view
colormapList = cellArray(colormapList);
viewSet(v,'colormaps',colormapList);

% display what we are doing
if verbose, disp(sprintf('(mlrAdjustGUI) Adding colormaps: %s',cellToCommaDelimited(colormapList))); end

%%%%%%%%%%%%%%%%%%%%%
%%%   addControl  %%%
%%%%%%%%%%%%%%%%%%%%%
function retval = addControl(f,args,uiControls,v,verbose)

% check length of arguments
if length(args) < 1
  disp(sprintf('(mlrAdjustGUI:addControl) Requires at least arguments: controlTag'));
  return
else
  % name the arguments
  controlTag = args{1};
  controlProperties = {args{2:end}};
  if isodd(length(controlProperties) )
    disp(sprintf('(mlrAdjustGUI:addControl) Properties must all have a matching property value'));
    return
  end
end

% check for panel
panelHandle = [];
for iArg = 1:2:length(controlProperties)
  if strcmp(lower(controlProperties{iArg}),'panel')
    % get panel name
    panelName = controlProperties{iArg+1};
    % ge tpanel handle
    panelHandle = viewGet(v,'panelHandle',panelName);
    if isempty(panelHandle)
      disp(sprintf('(mlrAdjustGUI) Could not find panel: %s',panelName));
    end
    % remove panel from controlProperties
    controlProperties = {controlProperties{1:iArg-1} controlProperties{iArg+2:end}};
    break;
  end
end

% check to see if it has already been added
if ~isempty(getHandle(controlTag,uiControls))
  disp(sprintf('(mlrAdjustGUI:addControl) Already added menu item: %s',controlTag));
  return
end

% get the gui data
h = guidata(f);

if ~isempty(panelHandle)
  % add to panel
  h.(controlTag)=uicontrol(f,'Parent',panelHandle);
else
  % add to main GUI
  h.(controlTag)=uicontrol(f);
end
set(h.(controlTag),'unit','normalized');
% add all the properties
for i = 1:2:length(controlProperties)
  set(h.(controlTag),controlProperties{i},controlProperties{i+1});
end

% store the guidata
guidata(f,h);

retval=h.(controlTag);

%%%%%%%%%%%%%%%%%%%
%%%   addPanel  %%%
%%%%%%%%%%%%%%%%%%%
function addPanel(f,args,v,verbose)

% first argument is panel name
if length(args)>1
  panelName = args{1};
else
  % default panel name
  panelName = 'Default panel';
end

% second argument should be percent size
% which is the height (in percent of the
% full figure height) - panels will get stuck
% on to the right hand side of the figure.
if length(args)>=2
  percentSize = args{2};
  if ~isnumeric(percentSize) || (length(percentSize)~=1) || (percentSize<0) || (percentSize>1)
    disp(sprintf('(mlrAdjustGUI) Panel must specify a percent size between 0 and 1'));
    return
  end
else
  % default to 25% percentSize
  percentSize = .25;
end

% check if this is the first panel, just so we can decide
% whether to put a separator above the menu item
mrGlobals;
if ~isfield(MLR,'panels') || (length(MLR.panels) == 0)
  separator = 'on';
else
  separator = 'off';
end

% add an invisible panel of the appropriate size to the GUI
mlrGuiSet(v,'addPanel',panelName,percentSize);
% add a menu item that allows you to turn on and off the panel
mlrAdjustGUI(v,'add','menu',panelName,'/View/Remove All Overlays','Callback',@mlrAdjustGUIPanel,'Checked','on','separator',separator);

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    mlrAdjustGUIPanel    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mlrAdjustGUIPanel(hObject,eventdata)

% get view
v = viewGet(getfield(guidata(hObject),'viewNum'),'view');

% get panelName
panelName = get(hObject,'Label');

% Check if we are on or off
if strcmp(get(hObject,'Checked'),'on')
  % hide panel
  mlrGuiSet(v,'hidePanel',panelName);
  % turn off check
  mlrAdjustGUI(v,'set',panelName,'Checked','off');
else
  % dispaly panel
  mlrGuiSet(v,'dispPanel',panelName);
  % turn off check
  mlrAdjustGUI(v,'set',panelName,'Checked','on');
end

%%%%%%%%%%%%%%%%%%%%%
%%%   addAxes  %%%
%%%%%%%%%%%%%%%%%%%%%
function addAxes(f,args,plotAxes,verbose)

% check length of arguments
if length(args) < 1
  disp(sprintf('(mlrAdjustGUI:addAxes) Requires at least arguments: axesTag'));
  return
else
  % name the arguments
  axesTag = args{1};
  axesProperties = {args{2:end}};
  if isodd(length(axesProperties) )
    disp(sprintf('(mlrAdjustGUI:addAxes) Properties must all have a matching property value'));
    return
  end
end

% check to see if it has already been added
if ~isempty(getHandle(axesTag,plotAxes))
  disp(sprintf('(mlrAdjustGUI:addAxes) Already added menu item: %s',axesTag));
  return
end

% get the gui data
h = guidata(f);

h.(axesTag)=axes('parent',f);
% add all the properties
for i = 1:2:length(axesProperties)
  set(h.(axesTag),axesProperties{i},axesProperties{i+1});
end

guidata(f,h);


%%%%%%%%%%%%%%%%%%
%%%   placeMenu  %%%
%%%%%%%%%%%%%%%%%%
function placeMenu(args,menuControls,mode,verbose,hFigure)

% check length of arguments
if length(args) < 2
  disp(sprintf('(mlrAdjustGUI:placeMenu) Requires 2 arguments: menuName, menuLocation'));
  return
else
  % name the arguments
  menuName = args{1};
  menuLocation = args{2};
  menuProperties = {args{3:end}};
  if isodd(length(menuProperties) )
    disp(sprintf('(mlrAdjustGUI:placeMenu) Properties must all have a matching property value'));
    return
  end
end

% go look for the item location
h = getHandle(menuLocation,menuControls);
% if not found, then print warning, return
if isempty(h)
  disp(sprintf('(mlrAdjustGUI:placeMenu) Could not find menu location: %s',menuLocation));
  return
end

% check to see if it has already been added
if strcmp(mode,'add') && ~isempty(getHandle(menuName,menuControls))
  disp(sprintf('(mlrAdjustGUI:placeMenu) Already added menu item: %s',menuName));
  return
end

% check to see if the location has a / on the end of it, which
% means to add it *underneath* the location specified
if menuLocation(end) == '/'
  switch(mode)
    case 'add'
      % add the menu to the parent
      hAdded = uimenu(h,'Label',menuName);
    case 'move'
      %change the menu's parent
      hAdded = menuName;
      set(hAdded,'Parent',h);
  end
      
  % and reorder to top
  hChildren = get(h,'Children');
  hChildren = [hChildren(hChildren~=hAdded); hAdded];
  set(h,'Children',hChildren);
else
  
  % get the parent
  hParent = get(h,'Parent');

  switch(mode)
    case 'add'
      % add the menu to the parent
      hAdded = uimenu(hParent,'Label',menuName);
    case 'move'
      %change the menu's parent
      hAdded = menuName;
      set(hAdded,'Parent',hParent);
  end

  % reorder the children so that the item created is below the menuLocation
  hChildren = get(hParent,'Children');
  % remove hAdded 
  hChildren = hChildren(hChildren~=hAdded);
  hChildren = [hChildren(1:find(hChildren==h)-1);hAdded;hChildren(find(hChildren==h):end)];
  % now set the children so that everything will be in the right order
  set(hParent,'Children',hChildren');
 
end

% add all the properties
for i = 1:2:length(menuProperties)
  set(hAdded,menuProperties{i},menuProperties{i+1});
  %if there is a tag property, add the new menu handle to the figure data
  if strcmp(mode,'add') && strcmp(menuProperties{i},'tag')
    handles = guidata(hFigure);
    if ~ismember(menuProperties{i+1},fieldnames(handles))
      handles.(menuProperties{i+1})=hAdded;
      guidata(hFigure,handles)
    end
  end
end

% display what we have done
if verbose
  switch(mode)
    case 'add'
      disp(sprintf('(mlrAdjustGUI:placeMenu) Added menu: %s',get(hAdded,'label')));
    case 'move'
      disp(sprintf('(mlrAdjustGUI:placeMenu) Moved menu: %s',get(hAdded,'label')));
  end
end

%%%%%%%%%%%%%%%%%%%%%
%%%   removeMenu  %%%
%%%%%%%%%%%%%%%%%%%%%
function removeMenu(args,menuControls,verbose)

% check length of arguments
if length(args) < 1
  disp(sprintf('(mlrAdjustGUI:removeMenu) Requires 1 argument: menuLocation'));
  return
else
  % name the arguments
  menuLocation = args{1};
end

% go look for the item location
h = getHandle(menuLocation,menuControls);
% if not found, then print warning, return
if isempty(h)
  disp(sprintf('(mlrAdjustGUI:removeMenu) Could not find menu location: %s',menuLocation));
  return
end

set(h,'visible','off')
%delete(h);

% display what we have done
if verbose, disp(sprintf('(mlrAdjustGUI:removeMenu) Removed menu: %s',menuLocation)); end

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  setItemProperty   %%
%%%%%%%%%%%%%%%%%%%%%%%%%
function setItemProperty(args,uiControls,menuControls,plotAxes,verbose)

% check arguments
if length(args) ~= 3
  disp(sprintf('(mlrAdjustGUI:setItemProperty) Requires 3 arguments: controlName, propertyName, propertyValue'));
  return
else
  % name the arguments
  controlName = args{1};
  propertyName = args{2};
  propertyValue = args{3};
end

% go look for the control
h = getHandle(controlName,uiControls,menuControls,plotAxes);

% if not found, then print warning, return
if isempty(h),
  disp(sprintf('(mlrAdjustGUI:setItemProperty) Could not find control: %s',controlName));
  return
end

controlType = get(h,'type');
if ~isempty(propertyName)
    propertyName = lower(propertyName);
end    

if strcmp(propertyName,'location')
  if ~strcmp(get(h,'type'),'uimenu')
    disp(['(mlrAdjustGUI) Cannot change location property for ' controlName]);
    return;
  end
  placeMenu({h,propertyValue},menuControls,'move',verbose);
  
else %if the property is not 'location', then we just set the property using set

  % check if the property exists
  fieldNames = lower(fieldnames(get(h)));

  if ~ismember(propertyName,fieldNames)
    % if not, warn and continue
    if isstr(propertyName)
        disp(sprintf('(mlrAdjustGUI) *** Could not find property %s of %s %s ***',propertyName,controlType,controlName));
    else
      disp(sprintf('(mlrAdjustGUI) *** Could not find correct property of %s %s. Did you pass in a string indicating a property to set ***',controlType,controlName));
    end
    return
  end

  % if it does, make sure we have a valid property to set it to
  %if ~any(strcmp(propertyName,{'Callback','TooltipString'}))
  validValues = set(h,propertyName);
  if ~isempty(validValues)
    tf = false;
    for j = 1:length(validValues) 
      if isequal(propertyValue,validValues{j}),tf = true;end
    end
    if ~tf
      disp(sprintf('(mlrAdjustGUI) Property value:'))
      disp(propertyValue);
      disp(sprintf('is invalid for %s %s property %s. Valid values are: ',controlType,controlName,propertyName));
      disp(validValues);
      return
    end
  end

  % if we got here, then the value is ok so set it
  if verbose
    disp(sprintf('(mlrAdjustGUI) Setting property %s of %s %s',propertyName,controlType,controlName));
  end;
  set(h,propertyName,propertyValue);
end
%%%%%%%%%%%%%%%%%%%%%%%
%%%  getUiControls  %%%
%%%%%%%%%%%%%%%%%%%%%%%
function uiControls = getUiControls(f,verbose)

if nargin == 1,verbose = 0;end
uiControls = [];

% get the gui data
h = guidata(f);
itemNames = fieldnames(h);

% for each guidata, check if it is a handle and
% not a menu item
for i = 1:length(itemNames)
  for j = 1:length(h.(itemNames{i}))
    if ishandle(h.(itemNames{i})(j))
      % check if it is not a menu item
      if isequal(get(h.(itemNames{i})(j),'Type'),'uicontrol')
	% keep track of this one.
	uiControls(end+1).tag = itemNames{i};
	uiControls(end).fieldNum = j;
	uiControls(end).h = h.(itemNames{i})(j);
      end
    end
  end
end

% print out info
if verbose
  for i = 1:length(uiControls)
    if uiControls(i).fieldNum == 1
      disp(sprintf('(mlrAdjustGUI) Found uicontrol: %i:%s',i,uiControls(i).tag));
      disp(get(uiControls(i).h,'Callback'));
    else
      disp(sprintf('(mlrAdjustGUI) Found uicontrol: %i:%s %i',i,uiControls(i).tag,uiControls(i).fieldNum));
    end
  end
end

%%%%%%%%%%%%%%%%%%%%%%%
%%%  getAxes  %%%
%%%%%%%%%%%%%%%%%%%%%%%
function plotAxes = getAxes(f,verbose)

if nargin == 1,verbose = 0;end
plotAxes = [];

% get the gui data
h = guidata(f);
itemNames = fieldnames(h);

% for each guidata, check if it is a handle and
% not a menu item
for i = 1:length(itemNames)
  for j = 1:length(h.(itemNames{i}))
    if ishandle(h.(itemNames{i})(j))
      % check if it is not a menu item
      if isequal(get(h.(itemNames{i})(j),'Type'),'axes')
	% keep track of this one.
	plotAxes(end+1).tag = itemNames{i};
	plotAxes(end).fieldNum = j;
	plotAxes(end).h = h.(itemNames{i})(j);
      end
    end
  end
end

% print out info
if verbose
  for i = 1:length(plotAxes)
    if plotAxes(i).fieldNum == 1
      disp(sprintf('(mlrAdjustGUI) Found axes: %i:%s',i,plotAxes(i).tag));
      disp(get(uiControls(i).h,'Callback'));
    else
      disp(sprintf('(mlrAdjustGUI) Found axes: %i:%s %i',i,plotAxes(i).tag,plotAxes(i).fieldNum));
    end
  end
end


%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  getMenuControls  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%
function menuControls = getMenuControls(f,menuControls)

% function can be called recursively, if being called
% at top level, default lists.
if nargin == 1
  menuControls = {};
  topLabel = [];
else
  topLabel = menuControls(end).fullLabel;
end

% get list of all children
c = get(f,'Children');
% see if there are anu menus
for i = 1:length(c)
  % if it is a menu item
  if isequal(get(c(i),'Type'),'uimenu')
    % then get its label
    menuControls(end+1).fullLabel = sprintf('%s/%s',topLabel,get(c(i),'Label'));
    menuControls(end).label = get(c(i),'Label');
    menuControls(end).tag = get(c(i),'Tag');
    menuControls(end).h = c(i);
    % and go look for submenus
    menuControls = getMenuControls(c(i),menuControls);
  end
end

