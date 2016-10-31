% mylegend.m
%
%      usage: mylegend({names},{symbols},pos,<a=axisNum>,<Location=NorthEast>)
%         by: justin gardner
%       date: 10/19/04
%    purpose: legend that puts up correct symbols too
%
%             pos is an optional argument that specifies
%             position of the legend--
%             1=top right 2=top left 3=bottom left 4=bottom right
%             Location can be specified using the builtin legend conventions
%             mylegend({'name 1','name 2'},{'ko','r+'},'Location=South');
%
%             You can specify the axis to make the legend for with the a arg
%             mylegend({'name 1','name 2'},{'ko','r+'},'a',gca);
%
%       e.g.: mylegend({'name 1','name 2'},{'ko','r+'})
%             mylegend({'name 1','name 2'},2);
%             mylegend({'name 1','name 2'});
%             you can also specify custom colors
%             mylegend({'name 1','name 2'},{{'ko' [0.8 0.8 0.8]},{'r+' [0.7 0 0]}});
%             or specify custom symbol types
%             mylegend({'name 1','name 2'},{{'ko' 'MarkerFaceColor=k' 'MarkerSize=12'},{'rs' 'MarkerFaceColor=r' 'MarkerEdgeColor=k'}});
%
function hLegend = mylegend(varargin)

hLegend = [];
if nargin < 1,help mylegend,return,end

% parse arguments
[names symbols pos optionalArgs] = parseArgs(varargin);

% parse optionalArgs
a = [];
getArgs(optionalArgs,{'a=[]'});

% use gca if no axis specified
if isempty(a),a = gca;end

% plot all the symbols invisibly to get handles for them
h = plotSymbols(a,symbols);

% call legend function
if ~isempty(h)
  hLegend = legend(a,h,names,'Location',pos);
else
  hLegend = legend(a,names,'Location',pos);
end  

% now fixup the legend (set color of text appropriately and turn off tex interpreting)
fixLegend(hLegend,h,names,symbols);

%%%%%%%%%%%%%%%%%%%
%%   fixLegend   %%
%%%%%%%%%%%%%%%%%%%
function fixLegend(hLegend,h,names,symbols)

% get the children of the legend
hLegendChildren = get(hLegend,'Children');

% go through the children and set each text item
lastColor = [];
for i = 1:length(hLegendChildren)
  if strcmp(get(hLegendChildren(i),'Type'),'text')
    % turn off tex interpreter
    set(hLegendChildren(i),'Interpreter','none');
    % set the font
    set(hLegendChildren(i),'FontName','Helvetica');
    set(hLegendChildren(i),'FontAngle','oblique');
    % if we have plotted symbols, figure out color from matching symbol
    c = [];
    if ~isempty(h)
      % check which name is corresponds to
      whichSymbol = first(find(strcmp(get(hLegendChildren(i),'String'),names)));
      if ~isempty(whichSymbol)
	% get the color of the symbol
	c = get(h(whichSymbol),'Color');
      end
    % otherwise use last color
    else
      c = lastColor;
    end
    if ~isempty(c)
      % and set the text to this color
      set(hLegendChildren(i),'Color',c);
    end
  end
  % get color of last symbol
  lastColor = get(hLegendChildren(i),'Color');
end

%%%%%%%%%%%%%%%%%%%%%
%%   plotSymbols   %%
%%%%%%%%%%%%%%%%%%%%%
function h = plotSymbols(a,symbols)

h = [];
hold on
x=0;y=0;
for iSymbol = 1:length(symbols)
  % get this Symbol
  thisSymbol = cellArray(symbols{iSymbol});
  % get symbol from old style
  if (length(thisSymbol) >= 1) && (length(thisSymbol{1}) <= 4)
    Symbol = thisSymbol{1};
    thisSymbol = {thisSymbol{2:end}};
  end
  % get color from old style
  oldColor = [];
  if (length(thisSymbol) >= 1) && (isnumeric(thisSymbol{1}))
    oldColor = thisSymbol{1};
    thisSymbol = {thisSymbol{2:end}};
  end
  % parse what the user has set for this symbol
  optionalArgs = {'Color','MarkerSize','MarkerFaceColor','MarkerEdgeColor'};
  for i = 1:length(optionalArgs),eval(sprintf('%s=[];',optionalArgs{i}));end
  getArgs(thisSymbol,optionalArgs);
  % override with oldColor if passed in
  if ~isempty(oldColor) Color = oldColor;end
  % now plot it
  h(iSymbol) = plot(a,x,y,Symbol,'Visible','off'); 
  for i = 1:length(optionalArgs)
    if ~isempty(eval(optionalArgs{i}))
      set(h(iSymbol),optionalArgs{i},eval(optionalArgs{i}));
    end
  end
end

%%%%%%%%%%%%%%%%%%%
%%   parseArgs   %%
%%%%%%%%%%%%%%%%%%%
function [names symbols pos optionalArgs] = parseArgs(args)

% get number of arguments
nargs = length(args);

% default to empty
symbols = {};
pos = [];
optionalArgs = {};

% now parse the arguments depending on how many were passed in
if nargs == 1
  names = args{1};
elseif nargs == 2
  names = args{1};
  if isnumeric(args{2})
    pos = args{2};
  elseif iscell(args{2})
    symbols = args{2};
  else
    optionalArgs = {args{2}};
  end
elseif nargs >= 3
  names = args{1};
  if iscell(args{2})
    symbols = args{2};
    if isnumeric(args{3})
      pos = args{3};
      optionalArgs = {args{4:end}};
    else
      optionalArgs = {args{3:end}};
    end
  else
    optionalArgs = {args{2:end}};
  end
end

% if passed in numbers then convert to strings
if isnumeric(names)
  for i = 1:length(names)
    namesStr{i} = mlrnum2str(names(i),'sigfigs=-1');
  end
  names = namesStr;
end
  
% old style calling, in which position was 1 of 4 locations
if pos == 1
  pos = 'NorthEast';
elseif pos == 2
  pos = 'NorthWest';
elseif pos == 3
  pos = 'SouthWest';
elseif pos == 4
  pos = 'SouthEast';
end

% see if user specified a "Location" argument, then pass that
% back instead of position
Location = [];
[argNames argValues optionalArgs] = getArgs(optionalArgs,{'Location=[]'});
if ~isempty(Location)
  pos = Location;
end

if isempty(pos),pos = 'NorthEast';end
 
if isstr(pos) && ~any(strcmp(pos,{'North','South','East','West','NorthEast','NorthWest','SouthEast','SouthWest','NorthOutside','SouthOutside','EastOutside','WestOutside','NorthEastOutside','NorthWestOutside','SouthEastOutside','SouthWestOutside','Best','BestOutside'}))
  disp(sprintf('(mylegend) Unknown legend location %s - changing to NorthEast',pos));
  pos = 'NorthEast';
end

% check symbols to see if it is a cell array of color
for iSymbol = 1:length(symbols)
  % check if they are all colors
  isListOfColors = 1;
  if ~isnumeric(symbols{iSymbol}) || (length(symbols{iSymbol}) ~= 3)
    isListOfColors = 0;
    break
  end
  % if so, then fix up into a symbol list pairing
  if isListOfColors
    for iSymbol = 1:length(symbols)
      newSymbols{iSymbol}{1} = 'ks';
      newSymbols{iSymbol}{2} = symbols{iSymbol};
      newSymbols{iSymbol}{3} = 'MarkerFaceColor';
      newSymbols{iSymbol}{4} = symbols{iSymbol};
      
    end      
    symbols = newSymbols;
  end
end