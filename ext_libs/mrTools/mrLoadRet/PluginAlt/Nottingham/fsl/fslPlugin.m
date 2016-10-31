% fslPlugin.m
%
%        $Id$ 
%      usage: DefaultPlugin(action,<thisView>)
%         by: julien besle
%       date: 12/13/10
%
function retval = fslPlugin(action,thisView)

% check arguments
if ~any(nargin == [1 2])
  help fslPlugin
  return
end

switch action
  % return a help string
  case {'help','h','?'}
    retval = 'Adds FSL functionalities: (1) adds items in menus ''Overlays'' and ''ROI'' to apply FSL FNIRT warp spline coefficients to overlays and/or ROIs. (2) adds item in "Analysis" menu to correct EPI distortions in timeseries using FSL FUGUE. (3) enables FLOBS and TFCE options in GLM v2 Plugin. You need to have FSL installed on your system and work under Linux or MacOSX to use this plugin.';
  
  case {'install','i'}
  % check for a valid view
  if (nargin ~= 2) || ~isview(thisView)
     disp(sprintf('(fslPlugin) Need a valid view to install plugin'));
  else
  % % %     %first check if there is already an Overlay Menu
  % % %     hFig = viewGet(thisView,'fignum');
  % % %     hFigChildren = get(hFig,'Children');
  % % %     overlaysMenuInstalled = 0;
  % % %     for i = 1:length(hFigChildren)
  % % %       % if it is a menu item
  % % %       if isequal(get(hFigChildren(i),'Type'),'uimenu')&& strcmp(get(hFigChildren(i),'Label'),'Overlays')
  % % %         overlaysMenuInstalled=1;
  % % %         break;
  % % %       end
  % % %     end
  % % %     if ~overlaysMenuInstalled
  % % %       %if not, install it
      mlrAdjustGUI(thisView,'add','menu','overlaysMenu','/Analysis','label','Overlays','tag','overlaysMenu');
  % %     end
    %install fslApplyWarp overlays menu Item
    mlrAdjustGUI(thisView,'add','menu','Apply FSL FNIRT non-linear warps to overlays','/Overlays/','callback',@applyWarpOverlaysCallback,'tag','applyFnirtOverlayMenuItem');
    %install fslApplyWarp ROIs menu Item
    mlrAdjustGUI(thisView,'add','menu','Apply FSL FNIRT non-linear warps to ROIs','/ROI/Combine','callback',@applyWarpROICallback,'tag','applyFnirtRoiMenuItem');
    %install FUGUE ROIs menu Item
    mlrAdjustGUI(thisView,'add','menu','Correct EPI distortions using FSL FUGUE','/Analysis/Motion Compensation','callback',@fugueTseriesCallback,'tag','applyFnirtRoiMenuItem');
    retval = true;
  end

  otherwise
    disp(sprintf('(fslPlugin) Unknown command %s',action));
end

%------------------------- applyWarpOverlaysCallback Function ------------------------------%
function applyWarpOverlaysCallback(hObject,dump)
thisView = viewGet(getfield(guidata(hObject),'viewNum'),'view');
fslApplyWarpOverlays(thisView);

%------------------------- applyWarpROICallback Function ------------------------------%
function applyWarpROICallback(hObject,dump)
thisView = viewGet(getfield(guidata(hObject),'viewNum'),'view');
fslApplyWarpROI(thisView);

%------------------------- applyWarpROICallback Function ------------------------------%
function fugueTseriesCallback(hObject,dump)
thisView = viewGet(getfield(guidata(hObject),'viewNum'),'view');
fslFugueTseries(thisView);
