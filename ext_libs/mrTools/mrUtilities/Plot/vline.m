% vline.m
%
%   usage: vline(hpos,<linetype>,<axis>)
%      by: justin gardner
%    date: 01/21/99
% purpose: draws a vertical line on the current axis
%          returns handles 
%          hpos can be an array of positions
function h = vline(hpos,linetype,a)

if nargin < 2,linetype = 'k:';end
if nargin < 3,a = gca;end

ax = axis(a);
miny = ax(3);maxy = ax(4);
if isequal(get(a,'YScale'),'log')
  miny = min(min(get(a,'Ytick')),miny);
end

h = [];
for i = 1:length(hpos)
  hold(a,'on');
  h(i) = plot(a,[hpos(i) hpos(i)],[miny maxy],linetype);
end
