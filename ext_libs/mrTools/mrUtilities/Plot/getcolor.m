% getcolor.m
%
%      usage: getcolor(colornum,symbol)
%         by: justin gardner
%       date: 10/07/04
%    purpose: return a different color for a different
%             number, optionally concatenate plot symbol
%       e.g.: getcolor(1)
%             getcolor(1,'o')
%
function retval = getcolor(colornum,symbol)

if ((nargin ~= 1) && (nargin ~= 2))
  help getcolor
  return
end

% get colors to use
[edgecolor colors] = getedgecolor;

if (nargin == 1)
  retval = colors(mod(colornum,length(colors))+1);
else
  retval = [colors(mod(colornum,length(colors))+1) symbol];
end