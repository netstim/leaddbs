% getSmoothColor.m
%
%        $Id: getSmoothColor.m,v 1.2 2008/07/16 21:38:18 justin Exp $ 
%      usage: getSmoothColor(colorNum,totalColors,<colorMap>,<skipRange>)
%         by: justin gardner
%       date: 07/08/08
%    purpose: returns a color that smoothly varies over totalColors
%             colorNum is a number from 1 to totalColors specifying the desired color
%             totalColors is the total number of colors that can be returned
%             colorMap is the colormap function to use (e.g. 'hot', 'cool', 'hsv', ...)
%                 it can also be an nx3 matrix containing a color map
%             if skipRange is specified, then you can prevent colors from going all
%             white, by specifying that you want to use only the bottom 80% of colors
%             for instance (i.e. by specifying skipRange to be 0.8); If you want
%             to use the top 80% of colors, make skipRange = -0.8;
%
%       e.g.: 
%
%for i = 1:10
%  plot(i,1,'ko','MarkerFaceColor',getSmoothColor(i,10),'MarkerSize',16);hold on
%  text(i,1,sprintf('%i',i),'HorizontalAlignment','Center','Color',getSmoothColor(11-i,10));
%end
%
%
function color = getSmoothColor(colorNum,totalColors,colorMap,skipRange)

% default return color
color = [0.5 0.5 0.5];

% check arguments
if ~any(nargin == [1 2 3 4])
  help getSmoothColor
  return
end

% default arguments
if ieNotDefined('totalColors'), totalColors = 256;end
if ieNotDefined('colorMap'), colorMap = 'gray';end
if ieNotDefined('skipRange')
  if strcmp(colorMap,'gray')
    skipRange = 0.8;
  else
    skipRange = 1;
  end
end


% get length of table necessary, this is modified by skipRange
% since we generate a longer than necessary table to skip the top
% most colors.
if skipRange > 0
  tableLen = ceil(totalColors*((1-skipRange)+1));
else
  tableLen = ceil(totalColors*((1+skipRange)+1));
end  

% if colorMap is a string, then it means to use that colormap function
if isstr(colorMap)
  % convert colorMap name into a function handle
  if ~any(strcmp(colorMap,{'hsv','gray','pink','cool','bone','copper','flag'}))
    if ~exist(colorMap,'file')
      disp(sprintf('(getSmoothColor) Unknown colormap function %s',colorMap));
      return
    end
  end
  % convert name of colormap function to a function handle
  % and generate a table of the correct length
  colorMapFun = str2func(colorMap);
  colors = colorMapFun(tableLen);
% colormap is a passed in matrix.
elseif isnumeric(colorMap)
  % check dimensions of color map, should be an nx3
  if ~size(colorMap,2) == 3
    disp(sprintf('(getSmoothColor) Passed in color map should be an nx3 array or a string'));
    return
  end

  % now resample if necessary
  if size(colorMap,1) ~= tableLen
    if size(colorMap,1) == 1
      colors = repmat(colorMap,tableLen,1);
    else
      [x y] = meshgrid(1:3,0:1/(size(colorMap,1)-1):1);
      [xi yi] = meshgrid(1:3,0:1/(tableLen-1):1);
      colors = interp2(x,y,colorMap,xi,yi);
    end
  else
    colors = colorMap;
  end
end  

% flip table appropriately for skipping bottom-most colors
if skipRange <= 0
  colors = colors(end-totalColors+1:end,:);
end

% select out the right color
if (colorNum >= 1) & (colorNum <= totalColors)
  color = colors(colorNum,:);
else
  % out of bounds. Warn and return gray
  disp(sprintf('(getSmoothColor) Color %i out of bounds [1 %i]',colorNum,totalColors));
end


