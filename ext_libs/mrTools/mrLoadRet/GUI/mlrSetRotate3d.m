% mlrSetRotate3d.m
%
%      usage: mlrSetRotate3d(v,onoff)
%         by: justin gardner
%       date: 07/11/09
%    purpose: Turn on/off mlr surface free rotation, v is the MLR
%             view. onoff is 1 or 0 to turn on or off surface viewing
% 
%
function retval = mlrSetRotate3d(v,onoff)

% check arguments
if ~any(nargin == [2])
  help mlrSetRotate3d
  return
end

% convert onoff strings to 0 or 1
if isstr(onoff)
  if strcmp(onoff,'on')
    onoff = 1;
  else
    onoff = 0;
  end
end

% get the figure number
f = viewGet(v,'fignum');
if isempty(f),
  disp(sprintf('(mlrSetRotate3d) No active figure for passed in view'));
  return
end

% turn on
if onoff
  h = rotate3d(f);
  % set the callback so that we can update the sliders
  set(h,'ActionPostCallback',@myPostCallback);
  set(h,'Enable','on');
  % have to keep the view as a global so that it is accessible from
  % the callback
  global gRotateView;
  gRotateView = v;
 else
  % turn off
  rotate3d(f,'off');
end

%%%%%%%%%%%%%%%%%%%%%%%%
%%   myPostCallback   %%
%%%%%%%%%%%%%%%%%%%%%%%%
function myPostCallback(obj,evd)

global gRotateView;

% get the most current version of the view.
v = viewGet([],'view',viewGet(gRotateView,'viewNum'));

% get what the axis was rotated to
newRotTilt = round(get(evd.Axes,'View'));

% get the tilt
tilt = newRotTilt(2);
if tilt < 0,tilt = -tilt;else,tilt = 360-tilt;end  
% and set the slider
v = viewSet(v,'tilt',tilt);

% get the rotation
rot = newRotTilt(1);
if rot < 0,rot = -rot;else,rot = 360-rot;end  
% and set the slider
v = viewSet(v,'rotate',rot);
