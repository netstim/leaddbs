% modelCoranalFromGlmPlugin.m
%
%        $Id$ 
%      usage: DefaultPlugin(action,<thisView>)
%         by: julien besle
%       date: 12/13/10
%    purpose: Plugin function for modelCoranalFromGlm
%
function retval = modelCoranalFromGlmPlugin(action,thisView)

% check arguments
if ~any(nargin == [1 2])
  help modelCoranalFromGlmPlugin
  return
end

switch action
 case {'install','i'}
  % check for a valid view
  if (nargin ~= 2) || ~isview(thisView)
     disp(sprintf('(modelCoranalFromGlmPlugin) Need a valid view to install plugin'));
  else
    mlrAdjustGUI(thisView,'add','menu','Model Correlation Analysis from GLM','/Analysis/GLM Analysis','Callback',@modelCoranalFromGlmCallback);
    retval = true;
   end
 % return a help string
 case {'help','h','?'}
   retval = 'Simulate a correlation analysis using regressor parameter estimate from a GLM analysis !!!';
 otherwise
   disp(sprintf('(modelCoranalFromGlmPlugin) Unknown command %s',action));
end

%------------------------- modelCoranalFromGlmCallback Function ------------------------------%
function modelCoranalFromGlmCallback(hObject,dump)
thisView = viewGet(getfield(guidata(hObject),'viewNum'),'view');
modelCoranalFromGlm(thisView);
