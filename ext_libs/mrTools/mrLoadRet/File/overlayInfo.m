% overlayInfo.m
%
%        $Id$
%      usage: overlayInfo(v)
%         by: justin gardner
%       date: 10/01/07
%    purpose: 
%
function retval = overlayInfo(v)

% check arguments
if ~any(nargin == [1])
  help overlayInfo
  return
end

% get the current overlay
o = viewGet(v,'overlay');

if isempty(o)
  disp(sprintf('(overlayInfo) No current overlay'));
  return
end

% get the fields and set some fields to not print
overlayFields = fieldnames(o);
ignoreFields = {'colormap','data','params'};

% set up the paramsInfo to display
paramsInfo = {};
for i = 1:length(overlayFields)
  if ~ismember(overlayFields{i},ignoreFields)
    paramsInfo{end+1} = {overlayFields{i},o.(overlayFields{i}),'editable=0'};
  end
end

% put in a button to put overlay into 
paramsInfo{end+1} = {'overlayData',[],'type=pushbutton','buttonString=Get overlay data','callback',@getOverlayData,'callbackArg',v,'Will set the variable overlayData in your current matlab workspace to the overlay.'};

mrParamsDialog(paramsInfo,'Overlay Info');

    
%%%%%%%%%%%%%%%%%%%%%%%%%%
% getOverlayData
%%%%%%%%%%%%%%%%%%%%%%%%%%
function retval = getOverlayData(v)

retval = [];
% get the current overlay
o = viewGet(v,'overlay');
disp(sprintf('(overlayInfo) Setting variable overlayData'));
assignin('base','overlayData',o.data);

