% color2RGB.m
%
%      usage: color = color2RGB(color)
%         by: justin gardner
%       date: 10/17/07
%    purpose: Converts a color name to an RGB value
%             With no arguments, gives a cell array
%             of color names
function color = color2RGB(color)

% check arguments
if ~any(nargin == [0 1])
  help color2RGB
  return
end

% get any color names that the user has set through mrSetPref - note
% that to avoid a calling loop we check whether we are being called by
% mrGetPref first before running this next section.
callingFunction = dbstack;
colorNames = {};
if (length(callingFunction)<1) || ~any(strcmp({callingFunction.name},'mrGetPref'))
  colorNames = mrGetPref('colorNames');
  % check to make sure that the colorNames are valid
  for iColorName = 1:length(colorNames)
    if ~iscell(colorNames{iColorName}) || (length(colorNames{iColorName}) ~= 2) || ~isstr(colorNames{iColorName}{1}) || ~isnumeric(colorNames{iColorName}{2}) || (length(colorNames{iColorName}{2})~=3)
      disp(sprintf('(color2RGB) colorNames preference is of an unrecognized fromat, so ignoring. colorNames should be a cell array of cell arrays, each element carrying a color name and color value triplet like: mrSetPref(''colorNames'',{{''mycolor'',[0.3 0.8 0.2]},{''myothercolor'',[0.9 0.1 0.2]}});'));
      colorNames = {};
      break;
    else
      % make sure that color names is a parsable name
      colorNames{iColorName}{1} = fixBadChars(colorNames{iColorName}{1});
    end
  end
end


% just return color names
if nargin == 0
  color = {'yellow','magenta','cyan','red','green','blue', ...
	   'orange','purple','white','black'};
  % add user defined colors
  for iColor = 1:length(colorNames)
    color{end+1} = colorNames{iColor}{1};
  end
  return
end

% convert color
if isstr(color)
  switch (color)
   case {'yellow','y'}, color = [1 1 0];
   case {'magenta','m'}, color = [1 0 1];
   case {'cyan','c'}, color = [0 1 1];
   case {'red','r'}, color = [1 0 0];
   case {'green','g'}, color = [0 1 0];
   case {'blue','b'}, color = [0 0 1];
   case {'orange','o'}, color = [255 165 0]/255;
   case {'purple','p'}, color = [160 32 240]/255;
   case {'white','w'}, color = [1 1 1];
   case {'black','k'}, color = [0 0 0];
   otherwise, 
    for iColor = 1:length(colorNames)
      % look for color in user supplied list
      if strcmp(lower(colorNames{iColor}{1}),lower(color));
	color = colorNames{iColor}{2};
	break;
      end
    end
  end % end switch statement
end

if isstr(color) || isempty(color),color = [1 1 1];end