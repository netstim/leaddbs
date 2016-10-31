% classificationAnalysisPlugin.m
%
%        $Id: classificationAnalysisPlugin.m 1969 2010-12-19 19:14:32Z julien $ 
%      usage: classificationAnalysisPlugin(action,<thisView>)
%         by: julien besle
%       date: 11/01/2011
%
function retval = classificationAnalysisPlugin(action,thisView)

% check arguments
if ~any(nargin == [1 2])
  help classificationAnalysisPlugin
  return
end

switch action
 case {'install','i'}
  % check for a valid view
  if (nargin ~= 2) || ~isview(thisView)
     disp(sprintf('(classificationAnalysisPlugin) Need a valid view to install plugin'));
  else
    %install menu Item
    mlrAdjustGUI(thisView,'add','menu','Classification','/Analysis/GLM analysis (v2)','callback',@classificationMenu_Callback);
    mlrAdjustGUI(thisView,'add','menu','Volume-Searchlight Classification','/Analysis/Classification/','callback',@searchlightClassification_Callback);
    mlrAdjustGUI(thisView,'add','menu','Surface-Searchlight Classification','/Analysis/Classification/','callback',@corticalSearchlightClassification_Callback);
    mlrAdjustGUI(thisView,'add','menu','ROI Classification','/Analysis/Classification/','callback',@roiClassification_Callback);
    mlrAdjustGUI(thisView,'add','menu','ROI GLM Classification','/Analysis/Classification/','callback',@glmClassification_Callback);

    retval = true;
   end
 % return a help string
 case {'help','h','?'}
   retval = 'Adds an item in Menu ''Analysis'' to do either spotlight or ROI based classification analysis';
 otherwise
   disp(sprintf('(classificationAnalysisPlugin) Unknown command %s',action));
end

% --------------------------------------------------------------------
function classificationMenu_Callback(hObject, eventdata)


% --------------------------------------------------------------------
function searchlightClassification_Callback(hObject, eventdata)

view = viewGet(getfield(guidata(hObject),'viewNum'),'view');
view = searchlightClassification(view);

function corticalSearchlightClassification_Callback(hObject, eventdata)

view = viewGet(getfield(guidata(hObject),'viewNum'),'view');
view = corticalSearchlight(view);

function roiClassification_Callback(hObject, eventdata)

view = viewGet(getfield(guidata(hObject),'viewNum'),'view');
view = roiClassification(view);

function glmClassification_Callback(hObject, eventdata)

view = viewGet(getfield(guidata(hObject),'viewNum'),'view');
view = eventRelatedROIClassification(view);