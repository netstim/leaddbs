% mrDispAsOverlay.m
%
%      usage: mrDispAsOverlay()
%         by: justin gardner
%       date: 04/04/07
%    purpose: 
%
function retval = mrDispAsOverlay()

% check arguments
if ~any(nargin == [1])
  help mrDispAsOverlay
  return
end

% create the parameters for the overlay
dateString = datestr(now);
quick.name = 'quick';
quick.function = '';
quick.groupName = 
quick.reconcileFunction = '';
quick.data = cell(1,viewGet(view,'nScans'));
quick.date = dateString;
quick.params = params;
quick.range = [0 1];
quick.clip = [0 1];
% colormap is made with a little bit less on the dark end
quick.colormap = hot(312);
quick.colormap = r2.colormap(end-255:end,:);
quick.alpha = 1;
quick.colormapType = 'setRangeToMax';
quick.interrogator = 'eventRelatedPlot';
