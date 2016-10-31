%
%        $Id:$ 
%      usage: fignum = mlrGetFignum(h);
%         by: justin gardner
%       date: 02/27/15
%    purpose: Get a figure number for a figure handle - this is to handle Mathworks
%             decision in R2014b to make figures handle objects and not numerics anymore...
%
function fignum = mlrGetFignum(h)

fignum = [];
% check arguments
if ~any(nargin == [1])
  help mlrGetFignum
  return
end

if isempty(h),return,end

% check if old type
if isnumeric(h)
  fignum = h;
else
  % check for number property
  if isprop(h,'Number') && ~isempty(h.Number)
    fignum = h.Number;
  else
    % can't figure out fignum
    disp(sprintf('(mlrGetFignum) Could not get a fignum for figure. Returning arbitray value'));
    fignum = round(rand*100);
  end
end

