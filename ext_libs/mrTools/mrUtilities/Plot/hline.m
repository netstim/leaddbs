% hline.m
%
%   usage: hline(vpos,<linetype>,<axis>)
%      by: justin gardner
%    date: 01/21/99
% purpose: draws a horizontal line on the current axis
%          returns handles 
%          vpos can be an array of positions
function h = hline(vpos,linetype,a)

if nargin < 2,linetype = 'k:';end
if nargin < 3,a = gca;end

ax = axis(a);
minx = ax(1);maxx = ax(2);
if isequal(get(a,'XScale'),'log')
  minx = min(get(a,'Xtick'));
end

h = [];
for i = 1:length(vpos)
  hold(a,'on');
  h(i) = plot(a,[minx maxx],[vpos(i) vpos(i)],linetype);
end

