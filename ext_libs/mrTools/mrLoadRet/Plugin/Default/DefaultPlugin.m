% DefaultPlugin.m
%
%        $Id:$ 
%      usage: DefaultPlugin(action,<v>)
%         by: justin gardner
%       date: 11/24/10
%    purpose: Plugin function for Test directory.
%
%             This follows a standard format in which
%             if action is set to 'install', then it adjusts the GUI using mlrAdjustGUI
%             to add the plugin functionaliy. If the action is set to 'help', then
%             it returns a string specifying help information for the plugin (usually
%             with a link to the wiki page which describes in more detail what functionality
%             the plugin offers).
%
%             Finally, this function should always be named with its parents direcory + the word Plugin.
%             e.g. if it is in the directory mrLoadRet/Plugin/Default, then it should be called DefaultPlugin.
%
%             That way the function mlrPlugin can find it appropriately.
%
function retval = DefaultPlugin(action,v)

% check arguments
if ~any(nargin == [1 2])
  help DefaultPlugin
  return
end

switch action
 case {'install','i'}
  % check for a valid view
  if (nargin ~= 2) || ~isview(v)
     disp(sprintf('(DefaultPlugin) Need a valid view to install plugin'));
  else
    % if the view is valid, then use mlrAdjustGUI to adjust the GUI for this plugin.
    
    % this installs a new menu item called 'Select Plugins' under /Edit/ROI with the
    % separator turned on above it. It sets the callback to selectPlugins defined below.
    mlrAdjustGUI(v,'add','menu','Select Plugins','/Edit/ROI','Callback',@selectPlugins,'Separator','on','tag','selectPluginMenuItem');

    % This is a command that could be used to install some default interrogators
    %mlrAdjustGUI(v,'add','interrogator',{'eventRelatedPlot','glmContrastPlot'});

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
   retval = 'This is an example plugin, it just installs a menu item to Select Plugins.';
 otherwise
   disp(sprintf('(DefaultPlugin) Unknown command %s',action));
end

%%%%%%%%%%%%%%%%%%%%%%%
%    selectPlugins    %
%%%%%%%%%%%%%%%%%%%%%%%
function selectPlugins(hObject,eventdata)

% This is an example callback function that gets called from the menu item "Select Plugins"
% It does not have to be a sub-function of the Plugin function as here, but could
% be a stand alone function

% code-snippet to get the view from the hObject variable. Not needed for this callback.
%v = viewGet(getfield(guidata(hObject),'viewNum'),'view');

% just run mlrPlugin to select which plugins to use
mlrPlugin;



