% pRFPlugin.m
%
%        $Id:$ 
%      usage: pRFPlugin(action,<v>)
%         by: justin gardner
%       date: 11/24/10
%    purpose: Plugin function for pRF directory.
%
function retval = pRFPlugin(action,v)

% check arguments
if ~any(nargin == [1 2])
  help pRFPlugin
  return
end

switch action
 case {'install','i'}
  % check for a valid view
  if (nargin ~= 2) || ~isview(v)
     disp(sprintf('(pRFPlugin) Need a valid view to install plugin'));
  else
    % if the view is valid, then use mlrAdjustGUI to adjust the GUI for this plugin.
    
    % this installs a new menu item called 'Select Plugins' under /Edit/ROI with the
    % separator turned on above it. It sets the callback to selectPlugins defined below.
    mlrAdjustGUI(v,'add','menu','pRF Analysis','/Analysis/Correlation Analysis','Callback',@callpRF,'Separator','off');

    % Install default interrogators
    mlrAdjustGUI(v,'add','interrogator',{'pRFFit'});

    % This is a command that could be used to install some default colormaps
    % that will show up when you do /Edit/Overlay
    %mlrAdjustGUI(v,'add','colormap','gray');

    % This is a command that could be used to set a property of an existing menu item
    %mlrAdjustGUI(v,'set','Plots/Mean Time Series','Separator','on');

    % return true to indicate successful plugin
    retval = true;
   end
 % return a help string
 case {'help','h','?'}
   retval = 'Runs population receptive field analysis.';
 otherwise
   disp(sprintf('pRFPlugin) Unknown command %s'));
end

%%%%%%%%%%%%%%%%%%%%%%%
%    selectPlugins    %
%%%%%%%%%%%%%%%%%%%%%%%
function callpRF(hObject,eventdata)

% code-snippet to get the view from the hObject variable. Not needed for this callback.
v = viewGet(getfield(guidata(hObject),'viewNum'),'view');
v = pRF(v);



