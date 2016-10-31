% viewGUIPlugin.m
%
%        $Id: viewGUIPlugin.m 1950 2010-12-18 10:12:48Z julien $ 
%      usage: viewGUIPlugin(action,<thisView>)
%         by: julien besle
%       date: 13/02/11
%    purpose: 
%

function retval = viewGUIPlugin(action,thisView)

% check arguments
if ~any(nargin == [1 2])
  help viewGUIPlugin
  return
end

switch action
 % return a help string
 case {'help','h','?'}
   retval = 'Reorganized menus, following the ''view'' logic';
   
 case {'install','i'}
  % check for a valid view
  if (nargin ~= 2) || ~isview(thisView)
     disp(sprintf('(viewGUIPlugin) Need a valid view to install plugin'));
  else
    
    %Create General Menu and move appropriate menu items
    mlrAdjustGUI(thisView,'add','menu','generalMenu','/File','label','General','tag','generalMenu');
    mlrAdjustGUI(thisView,'set','quitMenuItem','location','/General/');
    mlrAdjustGUI(thisView,'add','menu','quitWithoutSavingMenuItem','/General/','label','Quit without saving','tag','quitWithoutSaving','callback',@quitWithoutSavingMenuItem_Callback);
    mlrAdjustGUI(thisView,'add','menu','restartMenuItem','/General/','label','Restart','tag','restartMenuItem','callback',@restartMenuItem_Callback);
    mlrAdjustGUI(thisView,'add','menu','saveViewMenuItem','/General/','label','Save current view','tag','saveViewMenuItem','callback',@saveViewMenuItem_Callback);
    mlrAdjustGUI(thisView,'set','printMenuItem','location','/General/');
    mlrAdjustGUI(thisView,'set','exportImageMenuItem','location','/General/');
    mlrAdjustGUI(thisView,'set','selectPluginMenuItem','location','/General/');
    mlrAdjustGUI(thisView,'set','prefMenu','location','/General/');
    mlrAdjustGUI(thisView,'set','readmeMenuItem','location','/General/');
    mlrAdjustGUI(thisView,'set','editSessionMenuItem','location','/General/');
    mlrAdjustGUI(thisView,'set','saveSessionMenuItem','location','/General/');
    mlrAdjustGUI(thisView,'set','graphMenuItem','location','/General/');
    mlrAdjustGUI(thisView,'set','newWindowMenuItem','location','/General/');
    %rename menu items
    mlrAdjustGUI(thisView,'set','exportImageMenuItem','label','Export image(s)');
    %add/remove separators
    mlrAdjustGUI(thisView,'set','editSessionMenuItem','separator','on');
    mlrAdjustGUI(thisView,'set','readmeMenuItem','separator','off');
    mlrAdjustGUI(thisView,'set','printMenuItem','separator','off');
    mlrAdjustGUI(thisView,'set','exportImageMenuItem','separator','on');
    mlrAdjustGUI(thisView,'set','selectPluginMenuItem','separator','off');
    mlrAdjustGUI(thisView,'set','saveViewMenuItem','separator','on');
    mlrAdjustGUI(thisView,'set','quitMenuItem','separator','off');

    %Create Group Menu and move appropriate menu items
    mlrAdjustGUI(thisView,'add','menu','groupMenu','/General','label','Group','tag','groupMenu');
    mlrAdjustGUI(thisView,'set','deleteGroupMenuItem','location','/Group/');
    mlrAdjustGUI(thisView,'set','editGroupMenuItem','location','/Group/');
    mlrAdjustGUI(thisView,'set','infoGroupMenuItem','location','/Group/');
    mlrAdjustGUI(thisView,'set','importGroupMenuItem','location','/Group/');
    mlrAdjustGUI(thisView,'set','newGroupMenuItem','location','/Group/');
    %rename menu items
    mlrAdjustGUI(thisView,'set','deleteGroupMenuItem','label','Remove Selected');
    mlrAdjustGUI(thisView,'set','editGroupMenuItem','label','Edit');
    mlrAdjustGUI(thisView,'set','importGroupMenuItem','label','Import');
    mlrAdjustGUI(thisView,'set','newGroupMenuItem','label','New (empty)');
    
    %Create Scan Menu and move appropriate menu items
    mlrAdjustGUI(thisView,'add','menu','scanMenu','/Group','label','Scan','tag','scanMenu');
    mlrAdjustGUI(thisView,'set','scanViewInMlrVolMenuItem','location','/Scan/');
    mlrAdjustGUI(thisView,'set','editStimfileMenuItem','location','/Scan/');
    mlrAdjustGUI(thisView,'set','unlinkStimfileMenuItem','location','/Scan/');
    mlrAdjustGUI(thisView,'set','linkStimfileMenuItem','location','/Scan/');
    mlrAdjustGUI(thisView,'set','deleteScanMenuItem','location','/Scan/');
    mlrAdjustGUI(thisView,'set','pasteScanMenuItem','location','/Scan/');
    mlrAdjustGUI(thisView,'set','copyScanMenuItem','location','/Scan/');
    mlrAdjustGUI(thisView,'set','transformsMenu','location','/Scan/');
    mlrAdjustGUI(thisView,'set','dicomInfoMenuItem','location','/Scan/');
    mlrAdjustGUI(thisView,'set','infoScanMenuItem','location','/Scan/');
    mlrAdjustGUI(thisView,'set','exportTSeriesMenuItem','location','/Scan/');
    mlrAdjustGUI(thisView,'set','importTSeriesMenuItem','location','/Scan/');
    mlrAdjustGUI(thisView,'remove','menu''addScanMenuItem'); %this does the same thing as import
    %rename menu items
    mlrAdjustGUI(thisView,'set','deleteScanMenuItem','label','Remove...');
    mlrAdjustGUI(thisView,'set','addScanMenuItem','label','Add');
    mlrAdjustGUI(thisView,'set','importTSeriesMenuItem','label','Import');
    mlrAdjustGUI(thisView,'set','exportTSeriesMenuItem','label','Export');
    mlrAdjustGUI(thisView,'set','copyScanMenuItem','label','Copy...');
    mlrAdjustGUI(thisView,'set','pasteScanMenuItem','label','Paste');
    mlrAdjustGUI(thisView,'set','transformsMenu','label','Edit transform');
    %add/remove separators
    mlrAdjustGUI(thisView,'set','dicomInfoMenuItem','separator','off');
    mlrAdjustGUI(thisView,'set','infoScanMenuItem','separator','on');
   
    %Create Anatomy Menu and move appropriate menu items
    mlrAdjustGUI(thisView,'add','menu','anatomyMenu','/Scan','label','Anatomy','tag','anatomyMenu');
    mlrAdjustGUI(thisView,'set','Multiple base display','location','/Anatomy/');
    mlrAdjustGUI(thisView,'add','menu','removeAnatomyMenu','/Anatomy/','label','Remove','tag','removeAnatomyMenu');
    mlrAdjustGUI(thisView,'set','deleteAllBasesMenuItem','location','/Anatomy/Remove/');
    mlrAdjustGUI(thisView,'set','deleteManyBasesMenuItem','location','/Anatomy/Remove/');
    mlrAdjustGUI(thisView,'set','deleteBaseMenuItem','location','/Anatomy/Remove/');
    mlrAdjustGUI(thisView,'remove','menu','pasteBaseMenuItem');%will be replaced by 'duplicate Base'
    mlrAdjustGUI(thisView,'remove','menu','copyBaseMenuItem'); %will be replaced by 'duplicate Base'
    mlrAdjustGUI(thisView,'add','menu','duplicateBaseMenuItem','/Anatomy/','label','Duplicate selected','tag','duplicateBaseMenuItem','callback',@duplicateBaseMenuItem_Callback);
    mlrAdjustGUI(thisView,'set','baseTransformsMenu','location','/Anatomy/');
    mlrAdjustGUI(thisView,'set','editBaseMenuItem','location','/Anatomy/');
    mlrAdjustGUI(thisView,'set','infoBaseAnatomyMenuItem','location','/Anatomy/');
    mlrAdjustGUI(thisView,'remove','menu','exportAnatomyMenuItem'); %this is not implemented/would be the same as 'Save as ...
    mlrAdjustGUI(thisView,'set','Make plane','location','/Anatomy/');
    mlrAdjustGUI(thisView,'set','importSurfaceMenuItem','location','/Anatomy/');
    mlrAdjustGUI(thisView,'set','importFlatMenuItem','location','/Anatomy/');
    mlrAdjustGUI(thisView,'set','SaveAsMenuItem','location','/Anatomy/');
    mlrAdjustGUI(thisView,'set','saveAnatomyMenuItem','location','/Anatomy/');
    mlrAdjustGUI(thisView,'set','useCurrentScanMenuItem','location','/Anatomy/');
    mlrAdjustGUI(thisView,'set','loadFromVolumeMenuItem','location','/Anatomy/');
    mlrAdjustGUI(thisView,'set','loadAnatomyMenuItem','location','/Anatomy/');
    %rename menu items
    mlrAdjustGUI(thisView,'set','exportAnatomyMenuItem','label','Export');
    mlrAdjustGUI(thisView,'set','baseTransformsMenu','label','Edit transforms');
    mlrAdjustGUI(thisView,'set','deleteAllBasesMenuItem','label','All');
    mlrAdjustGUI(thisView,'set','deleteManyBasesMenuItem','label','Choose...');
    mlrAdjustGUI(thisView,'set','deleteBaseMenuItem','label','Selected');
    %add/remove separators
    mlrAdjustGUI(thisView,'set','importFlatMenuItem','separator','on');
    mlrAdjustGUI(thisView,'set','saveAnatomyMenuItem','separator','off');
    
    
    %move Analysis menu items
    mlrAdjustGUI(thisView,'add','menu','removeAnalysisMenu','/Analysis/','label','Remove','tag','removeAnalysisMenu');
    mlrAdjustGUI(thisView,'set','deleteAllAnalysisMenuItem','location','/Analysis/Remove/');
    mlrAdjustGUI(thisView,'set','deleteManyAnalysisMenuItem','location','/Analysis/Remove/');
    mlrAdjustGUI(thisView,'set','deleteAnalysisMenuItem','location','/Analysis/Remove/');
    mlrAdjustGUI(thisView,'set','pasteAnalysisMenuItem','location','/Analysis/');
    mlrAdjustGUI(thisView,'set','copyAnalysisMenuItem','location','/Analysis/');
    mlrAdjustGUI(thisView,'set','recomputeAnalysisMenuItem','location','/Analysis/');
    mlrAdjustGUI(thisView,'set','editAnalysisMenuItem','location','/Analysis/');
    mlrAdjustGUI(thisView,'set','EditAnalysisInfoMenuItem','location','/Analysis/');
    mlrAdjustGUI(thisView,'set','fileAnalysisMenu','location','/Analysis/');
    mlrAdjustGUI(thisView,'set','loadAnalysisMenuItem','location','/Analysis/');
    mlrAdjustGUI(thisView,'set','newAnalysisMenuItem','location','/Analysis/');
    %rename menu items
    mlrAdjustGUI(thisView,'set','newAnalysisMenuItem','label','New (empty)');
    mlrAdjustGUI(thisView,'set','copyAnalysisMenuItem','label','Copy');
    mlrAdjustGUI(thisView,'set','pasteAnalysisMenuItem','label','Paste');
    mlrAdjustGUI(thisView,'set','editAnalysisMenuItem','label','Edit');
    mlrAdjustGUI(thisView,'set','deleteAllAnalysisMenuItem','label','All');
    mlrAdjustGUI(thisView,'set','deleteManyAnalysisMenuItem','label','Choose...');
    mlrAdjustGUI(thisView,'set','deleteAnalysisMenuItem','label','Selected');
    mlrAdjustGUI(thisView,'set','fileAnalysisMenu','label','Save');
    mlrAdjustGUI(thisView,'set','saveAllAnalysisMenuItem','label','All');
    mlrAdjustGUI(thisView,'set','saveAnalysisMenuItem','label','Selected');
    mlrAdjustGUI(thisView,'set','recomputeAnalysisMenuItem','label','Edit Params & Recompute...');
    %add/remove separators
    mlrAdjustGUI(thisView,'set','recomputeAnalysisMenuItem','separator','off');
    mlrAdjustGUI(thisView,'set','motionCompMenu','separator','on');
    mlrAdjustGUI(thisView,'set','saveAnalysisMenuItem','separator','off');
    mlrAdjustGUI(thisView,'set','deleteAnalysisMenuItem','separator','off');


    %Create Overlays Menu and move appropriate menu items
    mlrAdjustGUI(thisView,'add','menu','overlaysMenu','/Analysis','label','Overlays','tag','overlaysMenu');
    mlrAdjustGUI(thisView,'set','viewMenu','location','/Overlays/');
    mlrAdjustGUI(thisView,'set','pasteOverlayMenuItem','location','/Overlays/');
    mlrAdjustGUI(thisView,'set','copyOverlayMenuItem','location','/Overlays/');
    mlrAdjustGUI(thisView,'set','editOverlayMenuItem','location','/Overlays/');
    mlrAdjustGUI(thisView,'set','overlayInfoMenuItem','location','/Overlays/');
    mlrAdjustGUI(thisView,'set','exportOverlayMenuItem','location','/Overlays/');
    mlrAdjustGUI(thisView,'set','importOverlayMenuItem','location','/Overlays/');
    mlrAdjustGUI(thisView,'set','fileOverlayMenu','location','/Overlays/');
    mlrAdjustGUI(thisView,'set','loadOverlayMenuItem','location','/Overlays/');
    %rename menu items
    mlrAdjustGUI(thisView,'set','exportOverlayMenuItem','label','Export');
    mlrAdjustGUI(thisView,'set','importOverlayMenuItem','label','Import');
    mlrAdjustGUI(thisView,'set','copyOverlayMenuItem','label','Copy...');
    mlrAdjustGUI(thisView,'set','pasteOverlayMenuItem','label','Paste');
    mlrAdjustGUI(thisView,'set','editOverlayMenuItem','label','Edit');
    mlrAdjustGUI(thisView,'set','viewMenu','label','Remove');
    mlrAdjustGUI(thisView,'set','deleteAllOverlayMenuItem','label','All');
    mlrAdjustGUI(thisView,'set','deleteManyOverlaysMenuItem','label','Choose');
    mlrAdjustGUI(thisView,'set','deleteOverlayMenuItem','label','Selected');
    mlrAdjustGUI(thisView,'set','fileOverlayMenu','label','Save');
    mlrAdjustGUI(thisView,'set','saveAllOverlayMenuItem','label','All');
    mlrAdjustGUI(thisView,'set','saveOverlayMenuItem','label','Selected');
    %add/remove separators
    mlrAdjustGUI(thisView,'set','deleteOverlayMenuItem','separator','off');
    mlrAdjustGUI(thisView,'set','saveOverlayMenuItem','separator','off');
    mlrAdjustGUI(thisView,'set','importOverlayMenuItem','separator','on');
    mlrAdjustGUI(thisView,'set','transformOverlaysMenuItem','separator','on');

    %move ROI menu items
    mlrAdjustGUI(thisView,'set','findCurrentROIMenuItem','location','/ROI/Show');
    mlrAdjustGUI(thisView,'set','findCurrentROIMenuItem','separator','off');
    mlrAdjustGUI(thisView,'set','deleteRoiMenu','location','/ROI/');
    mlrAdjustGUI(thisView,'set','convertRoiCoordsMenuItem','location','/ROI/');
    mlrAdjustGUI(thisView,'remove','menu','pasteRoiMenuItem'); %will be replaced by 'duplicate ROI'
    mlrAdjustGUI(thisView,'remove','menu','copyRoiMenuItem');  %will be replaced by 'duplicate ROI'
    mlrAdjustGUI(thisView,'add','menu','duplicateROIMenuItem','/ROI/','label','Duplicate selected','tag','duplicateROIMenuItem','callback',@duplicateRoiMenuItem_Callback);
    mlrAdjustGUI(thisView,'set','editRoiMenu','location','/ROI/');
    mlrAdjustGUI(thisView,'set','infoROIMenuItem','location','/ROI/');
    mlrAdjustGUI(thisView,'set','exportROIMenuItem','location','/ROI/');
    mlrAdjustGUI(thisView,'set','Import Freesurfer Label','location','/ROI/');
    mlrAdjustGUI(thisView,'set','Import Freesurfer Label','separator','off');
    mlrAdjustGUI(thisView,'set','importROIMenuItem','location','/ROI/');
    mlrAdjustGUI(thisView,'set','fileRoiMenu','location','/ROI/');
    mlrAdjustGUI(thisView,'set','loadFromVolumeDirectoryROIMenuItem','location','/ROI/');
    mlrAdjustGUI(thisView,'set','loadROIMenuItem','location','/ROI/');
    mlrAdjustGUI(thisView,'set','createRoiMenu','location','/ROI/');
    mlrAdjustGUI(thisView,'set','convertCorticalDepthRoiMenuItem','location','/ROI/Restrict');
    %rename menu items
    mlrAdjustGUI(thisView,'set','exportROIMenuItem','label','Export');
%     mlrAdjustGUI(thisView,'set','copyRoiMenuItem','label','Copy selected');
%     mlrAdjustGUI(thisView,'set','pasteRoiMenuItem','label','Paste');
    mlrAdjustGUI(thisView,'set','editRoiMenu','label','Edit');
    mlrAdjustGUI(thisView,'set','fileRoiMenu','label','Save');
    mlrAdjustGUI(thisView,'set','saveROIMenuItem','label','Selected');
    mlrAdjustGUI(thisView,'set','saveManyROIsMenuItem','label','Choose...');
    mlrAdjustGUI(thisView,'set','saveAllROIMenuItem','label','All');
    mlrAdjustGUI(thisView,'set','removeManyROIMenuItem','label','Choose...');
    mlrAdjustGUI(thisView,'set','convertRoiCoordsMenuItem','label','View/Convert Coordinates Space');
    mlrAdjustGUI(thisView,'set','convertCorticalDepthRoiMenuItem','label','Project through/Restrict to depth...');
    mlrAdjustGUI(thisView,'set','roiMenu','label','ROIs');
    mlrAdjustGUI(thisView,'set','addRoiMenu','label','Add Voxels');
    mlrAdjustGUI(thisView,'set','removeRoiMenu','label','Subtract Voxels');
    %add/remove separators
    mlrAdjustGUI(thisView,'set','saveROIMenuItem','separator','off');
    mlrAdjustGUI(thisView,'set','addRoiMenu','separator','on');
%     mlrAdjustGUI(thisView,'set','copyRoiMenuItem','separator','on');
    mlrAdjustGUI(thisView,'set','convertRoiMenuItem','separator','on');
    mlrAdjustGUI(thisView,'set','deleteRoiMenu','separator','off');    

    %remove unused menus
    mlrAdjustGUI(thisView,'remove','menu','fileMenu');
    mlrAdjustGUI(thisView,'remove','menu','editMenu');
    mlrAdjustGUI(thisView,'remove','menu','windowMenu');
    mlrAdjustGUI(thisView,'remove','menu','convertRoiMenuItem');
    
    % return 
    retval = true;
  end
  
  otherwise
   disp(sprintf('(viewGUIPlugin) Unknown command %s',action));
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Callbacks %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --------------------------------------------------------------------
function duplicateRoiMenuItem_Callback(hObject, eventdata)
handles = guidata(hObject);
viewNum = handles.viewNum;
rois = viewGet(viewNum,'roi');
for iRoi = 1:length(rois)
  view = viewSet(viewNum,'newROI',rois(iRoi));
end
% Select the last ROI
ROInum = viewGet(view,'numberofROIs');
if (ROInum > 0)
    view = viewSet(view,'currentROI',ROInum);
end
refreshMLRDisplay(viewNum);

% --------------------------------------------------------------------
function duplicateBaseMenuItem_Callback(hObject, eventdata)
handles = guidata(hObject);
viewNum = handles.viewNum;
viewSet(viewNum,'newBase',viewGet(viewNum,'baseAnatomy'));
refreshMLRDisplay(viewNum);
    
% --------------------------------------------------------------------
function restartMenuItem_Callback(hObject, eventdata)
mrGlobals;
handles = guidata(hObject);
viewNum = handles.viewNum;
v = MLR.views{viewNum};

mrQuit(1,v);
mrLoadRet;

% --------------------------------------------------------------------
function quitWithoutSavingMenuItem_Callback(hObject, eventdata)
mrGlobals;
handles = guidata(hObject);
viewNum = handles.viewNum;
v = MLR.views{viewNum};

mrQuit(0,v);

% --------------------------------------------------------------------
function saveViewMenuItem_Callback(hObject, eventdata)
mrGlobals;
handles = guidata(hObject);
viewNum = handles.viewNum;
v = MLR.views{viewNum};

mrSaveView(v);



