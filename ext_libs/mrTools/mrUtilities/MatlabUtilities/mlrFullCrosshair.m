% mlrFullCrosshair.m
%
%        $Id:$ 
%      usage: mlrFullCrosshair()
%         by: justin gardner
%       date: 03/12/15
%    purpose: An idiotic function required since matlab removed the full crosshair
%             mouse in Matlab 2015a (wtf Mathworks?) this returns fullcrosshair
%             or crosshair depending on version. Yuck
%
function retval = mlrFullCrosshair()

% check arguments
if ~any(nargin == [0])
  help mlrFullCrosshair
  return
end

if verLessThan('matlab','8.4')
  retval = 'fullcrosshair';
else
  retval = 'crosshair';
end

