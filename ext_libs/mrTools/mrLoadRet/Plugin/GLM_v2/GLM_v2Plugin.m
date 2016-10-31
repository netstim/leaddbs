% GLM_v2Plugin.m
%
%        $Id: GLM_v2Plugin.m 1950 2010-12-18 10:12:48Z julien $ 
%      usage: GLM_v2Plugin(action,<thisView>)
%         by: julien besle
%       date: 13/02/11
%    purpose: 
%

function retval = GLM_v2Plugin(action,thisView)

% check arguments
if ~any(nargin == [1 2])
  help GLM_v2Plugin
  return
end

switch action
  % return a help string
  case {'help','h','?'}
   retval =['GUI Improvements: - superimposition of several overlays\n',...
            '                  - average/projection over a range of cortical depths\n',...
            '                  - multiple ROI selection\n',... 
            '                  - folder mrLoadRet/Plugin/interrogatorFolder/interrogatorFunctions is scanned at MLR startup and functions to the default interrogator list',...
            '                  - folder mrLoadRet/Plugin/colormapFolder/colormapFunctions is scanned at MLR startup and adds functions to the default colormap list',...
            'Added functions: - New item in menu Edit/Scan menu to unlink stimfiles',...
            '                 - New item in menu ROIs to transform ROIs with pre-defined or custom transformation functions',...
            '                 - Improved GLM analysis that subsumes event-related and "GLM" analyses in a common GLM framework with statistical inference tests (see http://www.psychology.nottingham.ac.uk/staff/ds1/lab/doku.php?id=analysis:glm_statistics_in_mrtools',...
            '                 - New item in Plot that calls a 3D viewer',...
            '                 - New item in menu Overlays to combine/transform overlays using virtually any type of function',...
            '                 - Improved Edit/Overlay dialog',...
            '                 - Improved Export Images function',...
            '                 - New item in menu Analysis/Motioncomp to appy existing motion compensation parameters to another set of scans',...
            '                 - New item in menu edit/Base Anatomy/transforms to copy and paste sform',...
            ];
   retval ='Adds new functionalities to GUI, including improved GLM analysis';
   
  case {'install','i'}
  % check for a valid view
  if (nargin ~= 2) || ~isview(thisView)
     disp(sprintf('(GLM_v2Plugin) Need a valid view to install plugin'));
  else

    % new uicontrols and reposition old ones
    if viewGet(thisView,'baseType')>0
      corticalDepthVisibility = 'on';
    else
      corticalDepthVisibility = 'off';
    end    
    controlFontSize=10;
    checkFontSize=11;
    popupFontSize = 11;
    labelFontSize=12;
    sliderHeight=.025;
    textEditHeight=.035;
    checkBoxHeight=.025;
    boldTextLabelHeight=.03;
    
    %---------------------------- Group and Scan controls-----------------------------------------
    mlrAdjustGUI(thisView,'set','group','position',               [0.01    0.965   0.06   boldTextLabelHeight]);
    mlrAdjustGUI(thisView,'set','group','fontWeight','bold');
    mlrAdjustGUI(thisView,'set','groupPopup','position',          [0.08    0.955   0.2    0.035]);
    mlrAdjustGUI(thisView,'set','groupPopup','fontSize',popupFontSize);
    mlrAdjustGUI(thisView,'set','scan','position',                [0.01    0.925   0.05   boldTextLabelHeight]);
    mlrAdjustGUI(thisView,'set','scan','fontWeight','bold');
    mlrAdjustGUI(thisView,'set','scanSlider','position',          [0.07    0.92    0.16   sliderHeight ]);
    mlrAdjustGUI(thisView,'set','scanText','position',            [0.24    0.92    0.04   textEditHeight]);
    mlrAdjustGUI(thisView,'set','scanText','fontSize',controlFontSize);
    %---------------------------- Base controls  -------------------------------------------------
    mlrAdjustGUI(thisView,'set','baseImage','position',           [0.01    0.885   0.06   boldTextLabelHeight]);
    mlrAdjustGUI(thisView,'set','baseImage','fontWeight','bold');
    mlrAdjustGUI(thisView,'set','basePopup','position',           [0.07    0.875   0.21   0.035]);
    mlrAdjustGUI(thisView,'set','basePopup','fontSize',popupFontSize);
    mlrAdjustGUI(thisView,'set','slice','position',               [0.02    0.845   0.05   0.03 ]);
    mlrAdjustGUI(thisView,'set','slice','fontSize',checkFontSize);
    mlrAdjustGUI(thisView,'set','sliceSlider','position',         [0.08    0.845   0.15   sliderHeight ]);
    mlrAdjustGUI(thisView,'set','sliceText','position',           [0.24    0.845   0.04   textEditHeight]);
    mlrAdjustGUI(thisView,'set','sliceText','fontSize',controlFontSize);
    mlrAdjustGUI(thisView,'set','corticalDepth','position',       [0.02    0.815   0.07   0.055 ]);
    mlrAdjustGUI(thisView,'set','corticalDepth','fontSize',checkFontSize);
    mlrAdjustGUI(thisView,'set','corticalDepthSlider','position', [0.095   0.845   0.135  sliderHeight ]);
    mlrAdjustGUI(thisView,'set','corticalDepthSlider','SliderStep',min(1/viewGet(thisView,'corticalDepthBins')*[1 3],1));
    mlrAdjustGUI(thisView,'set','corticalDepthText','position',   [0.24    0.845   0.04   textEditHeight]);
    mlrAdjustGUI(thisView,'set','corticalDepthText','fontSize',controlFontSize);
    mlrAdjustGUI(thisView,'set','corticalDepthText','BackgroundColor', get(0,'defaultUicontrolBackgroundColor'));
    corticalDepthBins=viewGet(thisView,'corticaldepthbins');
    corticalDepth=round((corticalDepthBins-1)/2)/(corticalDepthBins-1);
    mlrAdjustGUI(thisView,'add','control','corticalMaxDepthSlider','SliderStep',min(1/(viewGet(thisView,'corticalDepthBins')-1)*[1 3],1),...
        'Callback',@corticalMaxDepthSlider_Callback,'visible',corticalDepthVisibility,'value',corticalDepth,...
                                     'style','slider','position', [0.095   0.82    0.135  sliderHeight ]);
    mlrAdjustGUI(thisView,'add','control','corticalMaxDepthText',...
        'Callback',@corticalMaxDepthText_Callback,'visible',corticalDepthVisibility,'String',num2str(corticalDepth),...
        'fontSize',controlFontSize,'BackgroundColor', get(0,'defaultUicontrolBackgroundColor'),...
                                        'style','edit','position',[0.24    0.815   0.04   textEditHeight ]);
    mlrGuiSet(thisView,'corticalmindepth',corticalDepth);
    mlrAdjustGUI(thisView,'add','control','linkMinMaxDepthCheck','style','checkbox','value',1,...
        'fontSize', checkFontSize,'visible',corticalDepthVisibility,...
        'String','Fix Depth Range','position',                    [0.095   0.79    0.13   checkBoxHeight]);
    mlrAdjustGUI(thisView,'set','sagittalRadioButton','position', [0.01    0.815   0.1    0.025]);
    mlrAdjustGUI(thisView,'set','sagittalRadioButton','fontSize',checkFontSize);
    mlrAdjustGUI(thisView,'set','axialRadioButton','position',    [0.11    0.815   0.07   0.025]);
    mlrAdjustGUI(thisView,'set','axialRadioButton','fontSize',checkFontSize);
    mlrAdjustGUI(thisView,'set','coronalRadioButton','position',  [0.19    0.815   0.1    0.025]);
    mlrAdjustGUI(thisView,'set','coronalRadioButton','fontSize',checkFontSize);

    % change multiAxis control position and fontSize
    mlrAdjustGUI(thisView,'set','axisSingle','position',  [0.01    0.788   0.1    checkBoxHeight]);
    mlrAdjustGUI(thisView,'set','axisSingle','fontSize',checkFontSize);
    mlrAdjustGUI(thisView,'set','axisMulti','position',  [0.11    0.788   0.1    checkBoxHeight]);
    mlrAdjustGUI(thisView,'set','axisMulti','fontSize',checkFontSize);
    mlrAdjustGUI(thisView,'set','axis3D','position',  [0.19    0.788   0.1    checkBoxHeight]);
    mlrAdjustGUI(thisView,'set','axis3D','fontSize',checkFontSize);

    mlrAdjustGUI(thisView,'set','baseGamma','string','Gamma');
    mlrAdjustGUI(thisView,'set','baseGamma','position',           [0.02    0.76    0.07   0.03 ]);
    mlrAdjustGUI(thisView,'set','baseGamma','fontSize',checkFontSize);
    mlrAdjustGUI(thisView,'set','baseGammaSlider','position',     [0.10    0.76    0.13   sliderHeight ]);
    mlrAdjustGUI(thisView,'set','baseGammaText','position',       [0.24    0.755   0.04   textEditHeight]);
    mlrAdjustGUI(thisView,'set','baseGammaText','fontSize',controlFontSize);
    mlrAdjustGUI(thisView,'set','rotate','position',              [0.02    0.73    0.06   0.03 ]);
    mlrAdjustGUI(thisView,'set','rotate','fontSize',checkFontSize);
    mlrAdjustGUI(thisView,'set','rotateSlider','position',        [0.09    0.73    0.14   sliderHeight ]);
    mlrAdjustGUI(thisView,'set','rotateText','position',          [0.24    0.725   0.04   textEditHeight]);
    mlrAdjustGUI(thisView,'set','rotateText','fontSize',controlFontSize);
    mlrAdjustGUI(thisView,'set','baseTilt','position',            [0.02    0.7     0.03   0.03 ]);
    mlrAdjustGUI(thisView,'set','baseTilt','fontSize',checkFontSize);
    mlrAdjustGUI(thisView,'set','baseTiltSlider','position',      [0.055   0.7     0.175  sliderHeight ]);
    mlrAdjustGUI(thisView,'set','baseTiltText','position',        [0.24    0.695   0.04   textEditHeight]);
    mlrAdjustGUI(thisView,'set','baseTiltText','fontSize',controlFontSize);

    %---------------------------- ROI controls  --------------------------------------------------
    mlrAdjustGUI(thisView,'set','ROI','position',                 [0.01    0.665   0.08    boldTextLabelHeight]);
    mlrAdjustGUI(thisView,'set','ROI','string','ROIs');
    mlrAdjustGUI(thisView,'set','ROI','fontWeight','bold');
    mlrAdjustGUI(thisView,'add','control','displayROILabels','style','checkbox','value',viewGet(thisView,'labelROIs'),...
        'callback',@displayROILabels_Callback,'fontSize', checkFontSize,'String',...
        'Display Labels','position',                              [0.09    0.665   0.18   checkBoxHeight]);
    switch(viewGet(thisView,'showROIs'))
      case {'selected'}
        roiDisplayMode = 1;
      case {'selected perimeter'}
        roiDisplayMode = 2;
      case {'group'}
        roiDisplayMode = 3;
      case {'group perimeter'}
        roiDisplayMode = 4;
      case {'all'}
        roiDisplayMode = 5;
      case {'all perimeter'}
        roiDisplayMode = 6;
      case 'hide'
        roiDisplayMode = 7;
    end
    mlrAdjustGUI(thisView,'add','control','roiDisplayMode','style','text','fontSize', checkFontSize,...
        'String','Display mode','position',                       [0.02    0.635   0.12   0.025]);
    mlrAdjustGUI(thisView,'add','control','roiDisplayModePopup','style','popupmenu','value',roiDisplayMode,...
        'callback',@roiDislayModePopup_Callback,'String',{'Selected ROIs (Voxels)' 'Selected ROIs (Perimeter)' 'ROIs group (Voxels)' 'ROIs group (Perimeter)' 'All ROIs (Voxels)' 'All ROIs (Perimeter)','Hide'},...
        'fontSize', popupFontSize,'position',                     [0.14    0.63    0.14   0.035]);
    mlrAdjustGUI(thisView,'set','roiPopup','position',            [0.02    0.51    0.26   0.125]);
    mlrAdjustGUI(thisView,'set','roiPopup','style','listbox');
    mlrAdjustGUI(thisView,'set','roiPopup','Max',2);  %allows multiselect
    mlrAdjustGUI(thisView,'set','roiPopup','callback',@roiList_Callback);
    mlrAdjustGUI(thisView,'set','roiPopup','value',viewGet(thisView,'currentRoi'));
    mlrAdjustGUI(thisView,'set','roiPopup','fontSize',controlFontSize);

    %---------------------------- Analysis and overlays controls  --------------------------------
    mlrAdjustGUI(thisView,'set','analysis','position',            [0.01    0.48    0.08   boldTextLabelHeight]);
    mlrAdjustGUI(thisView,'set','analysis','fontWeight','bold');
    mlrAdjustGUI(thisView,'set','analysisPopup','position',       [0.09    0.465   0.19   0.035]);
    mlrAdjustGUI(thisView,'set','analysisPopup','fontSize',popupFontSize);
    
    mlrAdjustGUI(thisView,'set','overlay','position',             [0.01    0.445   0.09   boldTextLabelHeight]);
    mlrAdjustGUI(thisView,'set','overlay','string','Overlays');
    mlrAdjustGUI(thisView,'set','overlay','fontWeight','bold');
    mlrAdjustGUI(thisView,'add','control','clipAcrossOverlays','style','checkbox','value',0,...
        'callback',@clipAcrossOverlays_Callback,'fontSize', checkFontSize,...
        'String','Clip across overlays','position',               [0.11    0.44    0.18   checkBoxHeight]);
    mlrAdjustGUI(thisView,'set','overlayPopup','position',        [0.02    0.275   0.26   0.16 ]);
    mlrAdjustGUI(thisView,'set','overlayPopup','style','listbox');
    mlrAdjustGUI(thisView,'set','overlayPopup','Max',2);  %allows multiselect
    mlrAdjustGUI(thisView,'set','overlayPopup','fontSize',controlFontSize);
    mlrAdjustGUI(thisView,'add','control','clippingOverlays','HorizontalAlignment','left','fontSize',labelFontSize,...
         'style','text','string','Clipping overlays','position',  [0.02    0.24    0.25   boldTextLabelHeight ]);
    mlrAdjustGUI(thisView,'add','control','clippingOverlaysListbox','callback',@clippingOverlaysCallback,...
        'fontSize',controlFontSize,'style','listbox','position',  [0.02    0.17    0.26   0.08 ]);
    mlrAdjustGUI(thisView,'set','overlayMin','position',          [0.02    0.135   0.04   0.03 ]);
    mlrAdjustGUI(thisView,'set','overlayMin','string','Min');
    mlrAdjustGUI(thisView,'set','overlayMin','fontSize',checkFontSize);
    mlrAdjustGUI(thisView,'set','overlayMinText','position',      [0.055   0.135   0.09   textEditHeight]);
    mlrAdjustGUI(thisView,'set','overlayMinText','fontSize',controlFontSize);
    mlrAdjustGUI(thisView,'set','mapMax','position',              [0.15    0.135   0.04   0.03 ]);
    mlrAdjustGUI(thisView,'set','mapMax','string','Max');
    mlrAdjustGUI(thisView,'set','mapMax','fontSize',checkFontSize);
    mlrAdjustGUI(thisView,'set','overlayMaxText','position',      [0.19    0.135   0.09   textEditHeight]);
    mlrAdjustGUI(thisView,'set','overlayMaxText','fontSize',controlFontSize);
    mlrAdjustGUI(thisView,'set','overlayMinSlider','position',    [0.02    0.11    0.26   sliderHeight ]);
    mlrAdjustGUI(thisView,'set','overlayMaxSlider','position',    [0.02    0.09    0.26   sliderHeight ]);
    mlrAdjustGUI(thisView,'set','alpha','position',               [0.02    0.055   0.05   0.03 ]);
    mlrAdjustGUI(thisView,'set','alpha','fontSize',checkFontSize);
    mlrAdjustGUI(thisView,'set','alphaSlider','position',         [0.08    0.055   0.14   sliderHeight ]);
    mlrAdjustGUI(thisView,'set','alphaText','position',           [0.23    0.055   0.05   textEditHeight]);
    mlrAdjustGUI(thisView,'set','alphaText','fontSize',controlFontSize);

    %---------------------------- Axes  --------------------------------
    mlrAdjustGUI(thisView,'set','colorbar','position',            [0.35    0.105   0.58   0.08 ]);
    mlrAdjustGUI(thisView,'set','colorbar','fontSize',10);
    mlrAdjustGUI(thisView,'add','axes','colorbarRightBorder',...
      'YaxisLocation','right','XTick',[],'box','off','color','none','position',  [0.929   0.105   0.001  0.08 ]);
    mlrAdjustGUI(thisView,'set','axis','position',                [0.285   0.195   0.71   0.8]);
    

    %----------------------------  Add menus --------------------------------

    % File menu
    mlrAdjustGUI(thisView,'set','exportImageMenuItem','callback',@exportImage_Callback);
    mlrAdjustGUI(thisView,'set','exportImageMenuItem','label','Export Images');

    % Edit menu
    mlrAdjustGUI(thisView,'add','menu','Unlink Stimfile','/Edit/Scan/Link Stimfile','callback',@unlinkStimfileMenuItem_Callback,'tag','unlinkStimfileMenuItem');
    mlrAdjustGUI(thisView,'set','/Edit/Scan/Link Stimfile','separator','on');
    mlrAdjustGUI(thisView,'add','menu','Copy sform','/Edit/Base Anatomy/Transforms/','callback',@copyBaseSformCallBack,'tag','copyBaseSformMenuItem');
    mlrAdjustGUI(thisView,'add','menu','Paste sform','/Edit/Base Anatomy/Transforms/','callback',@pasteBaseSformCallBack,'tag','pasteBaseSformMenuItem');
    mlrAdjustGUI(thisView,'set','editOverlayMenuItem','Callback',@editOverlayCallback);
    mlrAdjustGUI(thisView,'set','copyOverlayMenuItem','Callback',@copyOverlayCallback);
    mlrAdjustGUI(thisView,'set','copyOverlayMenuItem','label','Copy overlay(s)...');
    mlrAdjustGUI(thisView,'set','pasteOverlayMenuItem','Callback',@pasteOverlayCallback);
    mlrAdjustGUI(thisView,'set','copyScanMenuItem','Callback',@copyScanCallback);
    mlrAdjustGUI(thisView,'set','copyScanMenuItem','label','Copy scan(s)...');
    mlrAdjustGUI(thisView,'set','pasteScanMenuItem','Callback',@pasteScanCallback);

    % Analysis menu
    mlrAdjustGUI(thisView,'add','menu','Apply MotionComp Transforms','/Analysis/Motion Compensation/Slice Time Correction (only)','callback',@applyMotionCompTransformsCallBack,'tag','applyMotionCompTransformMenuItem');
    %mlrAdjustGUI(thisView,'remove','menu','eventRelatedMenuItem');
    mlrAdjustGUI(thisView,'set','glmMenuItem','Callback',@glmAnalysisCallback);
    mlrAdjustGUI(thisView,'set','glmMenuItem','label','GLM analysis (v2)');
    mlrAdjustGUI(thisView,'set','recomputeAnalysisMenuItem','Callback',@recomputeAnalysisCallback);

    % Overlay menu
    mlrAdjustGUI(thisView,'add','menu','overlaysMenu','/Analysis','label','Overlays','tag','overlaysMenu');
    %install menu Item
    mlrAdjustGUI(thisView,'add','menu','Combine/Transform Overlays','/Overlays/','callback',@combineOverlaysCallback,'tag','transformOverlaysMenuItem');


    %ROI menu
    %remove show group ROI menus and re-arrange ROI menus 
%     mlrAdjustGUI(thisView,'remove','menu','showGroupMenuItem');
%     mlrAdjustGUI(thisView,'remove','menu','showGroupPerimeterMenuItem');
%     mlrAdjustGUI(thisView,'remove','menu','setROIGroupMenuItem');
    mlrAdjustGUI(thisView,'set','undoRoiMenuItem','location','/ROI/Restrict');
    mlrAdjustGUI(thisView,'set','convertRoiMenuItem','separator','off');
    mlrAdjustGUI(thisView,'set','deleteRoiMenu','separator','on');
    mlrAdjustGUI(thisView,'set','restrictRoiMenuItem','label','Selected ROI(s)');
    mlrAdjustGUI(thisView,'set','deleteRoiMenuItem','label','Selected ROI(s)');
    mlrAdjustGUI(thisView,'set','saveROIMenuItem','label','Save selected ROI(s)');
    mlrAdjustGUI(thisView,'set','copyRoiMenuItem','label','Copy selected ROI(s)');
    mlrAdjustGUI(thisView,'set','pasteRoiMenuItem','label','Paste ROI(s)');
    mlrAdjustGUI(thisView,'set','editRoiMenuItem','label','Edit selected ROI(s)');
    %add functions
    mlrAdjustGUI(thisView,'add','menu','Single Voxels','/ROI/Create/Contiguous Voxels','callback',@createSingleVoxelsCallBack,'label','Single Voxels','tag','createSingleVoxelsRoiMenuItem','accelerator','T');
    mlrAdjustGUI(thisView,'add','menu','Single Voxels2','/ROI/Add/Contiguous Voxels','callback',@addSingleVoxelsCallBack,'label','Single Voxels','tag','addSingleVoxelsRoiMenuItem','accelerator','N');
    mlrAdjustGUI(thisView,'add','menu','Single Voxels3','/ROI/Subtract/Contiguous Voxels','callback',@removeSingleVoxelsCallBack,'label','Single Voxels','tag','removeSingleVoxelsRoiMenuItem','accelerator','U');
    mlrAdjustGUI(thisView,'add','menu','Transform','/ROI/Combine','callback',@transformROIsCallback,'label','Transform','tag','transformRoiMenuItem');

    %Plot menu
    %add 3D render viewer
    mlrAdjustGUI(thisView,'add','menu','3D Viewer','flatViewerMenuItem','callback',@renderIn3D,'tag','viewIn3DMenuItem');

    
    %---------------------------- Add colormaps and interrogators
  
    %get interrogators in the interrogatorFunctions directory
    interrogatorsFolder = [fileparts(which('GLM_v2Plugin')) '/interrogatorFunctions/'];
    interrogatorFiles =  dir([interrogatorsFolder '*.m']);
    if ~isempty(interrogatorFiles)
      interrogatorList = cell(1,length(interrogatorFiles));
      for iFile=1:length(interrogatorFiles)
         interrogatorList{iFile} = stripext(interrogatorFiles(iFile).name);
      end
      % install default interrogators
      % that will show up when you do /Edit/Overlay
      mlrAdjustGUI(thisView,'add','interrogator',interrogatorList);
    else
      disp('(interrogatorFolderPlugin) No interrogator function found in folder');
    end
    
    %get colormaps in the colormapFunctions directory
    colorMapsFolder = [fileparts(which('GLM_v2Plugin')) '/colormapFunctions/'];
    colorMapFiles =  dir([colorMapsFolder '*.m']);
    if ~isempty(colorMapFiles)
      colorMapList = cell(1,length(colorMapFiles));
      for iFile=1:length(colorMapFiles)
         colorMapList{iFile} = stripext(colorMapFiles(iFile).name);
      end
      % install default colormaps
      % that will show up when you do /Edit/Overlay
      mlrAdjustGUI(thisView,'add','colormap',colorMapList);
    else
      disp('(colormapFolderPlugin) No colormap function found in folder');
    end
   
    % return view
    retval = thisView;
  end

  otherwise
    disp(sprintf('(GLM_v2Plugin) Unknown command %s',action));
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Control Callbacks %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% -------------------------------------------------------------------- Cortical depth controls

% ----------------------------------- Executes on slider movement.
function corticalMaxDepthSlider_Callback(hObject, dump)

handles = guidata(hObject);
viewNum = handles.viewNum;
newMaxValue = get(hObject,'Value');
if get(handles.linkMinMaxDepthCheck,'value')
  minDepth = viewGet(viewNum,'corticalMinDepth');
  maxDepth = str2num(get(handles.corticalMaxDepthText,'string'));
  newMinValue = minDepth+newMaxValue-maxDepth;
  if newMinValue>=0 && newMinValue<=1 %if both values are in [0 1]
    viewSet(viewNum,'corticalMinDepth',newMinValue);
    viewSet(viewNum,'corticalMaxDepth',newMaxValue);
    drawnow;
    refreshMLRDisplay(viewNum);
  else
    set(hObject,'Value',maxDepth);
  end
else
  viewSet(viewNum,'corticalMaxDepth',newMaxValue);
  refreshMLRDisplay(viewNum);
end


% -------------------------------------------------------------------- 
function corticalMaxDepthText_Callback(hObject,dump)

handles = guidata(hObject);
viewNum = handles.viewNum;
newMaxValue = str2num(get(hObject,'String'));
if isempty(newMaxValue) %if the user just erased the value, get it from the slider and do nothing
  set(hObject,'String',num2str(get(handles.corticalMaxDepthSlider,'value')));
else %otherwise, set the new value in the view and the GUI
  if get(handles.linkMinMaxDepthCheck,'value')
    minDepth = viewGet(viewNum,'corticalMinDepth');
    maxDepth = get(handles.corticalMaxDepthSlider,'value');
    newMinValue = minDepth+newMaxValue-maxDepth;
    if newMinValue>=0 && newMinValue<=1 %if both values are in [0 1]
      viewSet(viewNum,'corticalMinDepth',newMinValue);
      viewSet(viewNum,'corticalMaxDepth',newMaxValue);
      refreshMLRDisplay(viewNum);
    else
      set(hObject,'string',num2str(maxDepth));
    end
  else
    viewSet(viewNum,'corticalMaxDepth',newMaxValue);
    refreshMLRDisplay(viewNum);
  end
end

% -------------------------------------------------------------------- Overlay controls
function clipAcrossOverlays_Callback(hObject,dump)
handles = guidata(hObject);
viewNum = handles.viewNum;

viewSet(viewNum,'clipAcrossOverlays',get(hObject,'value'));

%set overlay names in clipping box 
clippingOverlayList=viewGet(viewNum,'clippingOverlayList');
mlrGuiSet(viewNum,'clippingOverlays',unique(clippingOverlayList));

% set overlay min and max sliders
curClippingOverlay = viewGet(viewNum,'curClippingOverlay');
if ~isempty(curClippingOverlay)
  overlayClip = viewGet(viewNum,'overlayClip',curClippingOverlay);
  overlayRange = viewGet(viewNum,'overlayRange',curClippingOverlay);
  if ~isempty(overlayRange)
    mlrGuiSet(viewNum,'overlayMinRange',overlayRange);
    mlrGuiSet(viewNum,'overlayMaxRange',overlayRange);
  end
  if ~isempty(overlayClip)
    mlrGuiSet(viewNum,'overlayMin',overlayClip(1));
    mlrGuiSet(viewNum,'overlayMax',overlayClip(2));
  end
end

refreshMLRDisplay(viewNum);


% -------------------------------------------------------------------- Overlay controls
function clippingOverlaysCallback(hObject,dump)
handles = guidata(hObject);
viewNum = handles.viewNum;

clippingOverlay=viewGet(viewNum,'curClippingOverlay');
clip = viewGet(viewNum, 'overlayclip',clippingOverlay);
range = viewGet(viewNum, 'overlayrange',clippingOverlay);
mlrGuiSet(viewNum,'overlayMinRange',range);
mlrGuiSet(viewNum,'overlayMaxRange',range);
mlrGuiSet(viewNum,'overlayMin',clip(1));
mlrGuiSet(viewNum,'overlayMax',clip(2));

% -------------------------------------------------------------------- ROI controls
function roiList_Callback(hObject, dump)
handles = guidata(hObject);
viewNum = handles.viewNum;

roiNames = viewGet(viewNum,'roiNames');% I think we shouldn't have to get the names (but need to modify viewSet call)
% set the roi group
% viewSet(viewNum,'roiGroup',roiNames(get(hObject,'Value')));
viewSet(viewNum,'curROI',get(hObject,'Value'));
refreshMLRDisplay(viewNum);

function roiDislayModePopup_Callback(hObject, dump)
handles = guidata(hObject);
viewNum = handles.viewNum;
switch(get(hObject,'value'))
  case 1
    viewSet(viewNum,'showROIs','selected');
  case 2
    viewSet(viewNum,'showROIs','selected perimeter');
  case 3
    viewSet(viewNum,'showROIs','group');
  case 4
    viewSet(viewNum,'showROIs','group perimeter');
  case 5
    viewSet(viewNum,'showROIs','all');
  case 6
    viewSet(viewNum,'showROIs','all perimeter');
  case 7
    viewSet(viewNum,'showROIs','hide');
end
refreshMLRDisplay(viewNum);


% --------------------------------------------------------------------
function displayROILabels_Callback(hObject, dump)
handles = guidata(hObject);
viewNum = handles.viewNum;
viewSet(viewNum,'labelROIs',get(hObject,'Value'));
refreshMLRDisplay(viewNum);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% File Menu Callbacks %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --------------------------------------------------------------------
function exportImage_Callback(hObject, dump)
handles = guidata(hObject);
viewNum = handles.viewNum;
thisView = viewGet(viewNum,'view');

viewType = viewGet(thisView,'basetype');

if viewType<2
  nScans = viewGet(thisView,'nScans');
  paramsInfo = {};
  paramsInfo{end+1} = {'tifFileName', fullfile(viewGet(thisView,'homedir'),'Images/image.tif'),'Where to save the image montage file (TIF format)'};
  paramsInfo{end+1} = {'horizontalRange',[0 1],'minmax=[0 1]','Horizontal coordinates of the image cropped from the slice [left right] (normalized between 0 and 1)'};
  paramsInfo{end+1} = {'verticalRange',[0 1],'minmax=[0 1]','Vertical coordinates of the image cropped from the slice [top bottom] (normalized between 0 and 1)'};
  paramsInfo{end+1} = {'scanList', sprintf('[1:1:%d]',nScans), 'List of scans to export. Must be expressed in a string taht can be evaluated into a row vector'};

  if viewType
    paramsInfo{end+1} = {'depthMode', {'as displayed','one image per depth'}, '''As displayed'' = only one image is exported per scan, with the current depth settings, ''One image per depth''= exports one image per depth/scan within the displayed depth range'};
  else
    basedims = viewGet(thisView,'basedims');
    nSlices = basedims(viewGet(thisView,'basesliceindex'));

    paramsInfo{end+1} = {'sliceList', sprintf('[1:1:%d]',nSlices), 'List of slices to export. Must be expressed in a string that can be evaluated into a row vector'};
  end
  paramsInfo{end+1} = {'nRows', 1, 'minmax=[1 inf]','incdec=[-1 1]','Number of rows in the montage if multiple scans/slices/depths are exported'}; 
  paramsInfo{end+1} = {'exportImage', 0, 'type=pushbutton','','callback',@exportCallback,'callbackArg',thisView,'buttonString=Export Image','passParams=1','Previews and saves the montage without closing the menu'};

  % display dialog
  mrParamsDialog(paramsInfo,'Export Slices Parameters','modal=0');
else
    % if surface, would need to change the ouput of refreshMLRDisplay to get a 2D image rather than the coordinates of the surface...
    % but not sure whether some other parts of mrLoadRet use this output..
    mrWarnDlg('(GLM_v2Plugin) This function is not implemented for surfaces.');
end


% --------------------------------------------------------------------
function exportCallback(thisView,params)

%update the view
thisView = viewGet(thisView.viewNum,'view');

if viewGet(thisView,'basetype')
  switch(params.depthMode)
    case 'as displayed'
      sliceList = [];
    case 'one image per depth'
      depthBin = 1/(viewGet(thisView,'corticaldepthBins')-1);
      depth(1) = viewGet(thisView,'corticalmindepth');
      depth(2) = viewGet(thisView,'corticalmaxdepth');
      depth = sort(depth);
      sliceList = depth(1):depthBin:depth(2);
%       [~,sliceList] = ismember(round(1e5*(depth(1):depthBin:depth(2)))/1e5,round(1e5*(0:depthBin:1))/1e5);
  end
else
  sliceList = eval(params.sliceList);
end

scanList = eval(params.scanList);
mrSliceExport(thisView, [params.horizontalRange params.verticalRange], sliceList, params.tifFileName, params.nRows, scanList)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Edit Menu Callbacks %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --------------------------------------------------------------------
function copyOverlayCallback(hObject, dump)

mrGlobals;
thisView = viewGet(getfield(guidata(hObject),'viewNum'),'view');
clipboard = copyOverlay(thisView); %calls copyOverlay which asks which overlay and for which scans to copy an returns all the copied overlays
if ~isempty(clipboard)
  MLR.clipboard = clipboard;
end

% --------------------------------------------------------------------
function pasteOverlayCallback(hObject, dump)

mrGlobals;
thisView = viewGet(getfield(guidata(hObject),'viewNum'),'view');
if isfield(MLR,'clipboard') 
  pasteOverlay(thisView, MLR.clipboard);
  refreshMLRDisplay(thisView.viewNum);
else
  mrWarnDlg('(pasteOverlayCallback) There is no overlay to paste.')
end


% --------------------------------------------------------------------
function copyScanCallback(hObject, dump)

mrGlobals;
thisView = viewGet(getfield(guidata(hObject),'viewNum'),'view');
nScans = viewGet(thisView,'nScans');
if nScans>1
  scanList = selectInList(thisView,'scans','Select Scans to copy');
elseif nScans==1
  scanList = 1;
else
  mrWarnDlg('(copyScanCallback) There is no scan to copy');
  return;
end
scan=[];
cScan=0;
for iScan = scanList
  cScan = cScan+1;
  scan(cScan).scanParams = viewGet(thisView,'scanParams',iScan);
  scan(cScan).scanParams.fileName = [];
  % just save the filename instead of loading the whole tSeries
  % since some scans are very large
  %scan.tseries = loadTSeries(thisView,cScan,'all');
  scan(cScan).tseries = viewGet(thisView,'tseriesPathStr',iScan);
end
MLR.clipboard = scan;

% --------------------------------------------------------------------
function pasteScanCallback(hObject, dump)

mrGlobals;
thisView = viewGet(getfield(guidata(hObject),'viewNum'),'view');
if isfield(MLR,'clipboard') && isfield(MLR.clipboard,'tseries') && isfield(MLR.clipboard,'scanParams') 
  for iScan = 1:length(MLR.clipboard)
    if isscan(MLR.clipboard(iScan).scanParams)
    	thisView = saveNewTSeries(thisView,MLR.clipboard(iScan).tseries,MLR.clipboard(iScan).scanParams,MLR.clipboard(iScan).scanParams.niftiHdr);
    else
      mrWarnDlg(['(paste scan) Could not paste scan ' MLR.clipboard(iScan).tseries ' because its parameters are not valid'])
    end
  end
else
    mrWarnDlg('(paste scan) Cannot paste. Clipboard does not contain valid scans. Use Edit -> Scan -> Copy Scan.')
end

% -------------------------------------------------------------------- Unlink stim menu
function unlinkStimfileMenuItem_Callback(hObject, dump)
handles = guidata(hObject);
viewNum = handles.viewNum;
viewSet(viewNum, 'stimfilename', '');


% --------------------------------------------------------------------
function copyBaseSformCallBack(hObject, dump)

mrGlobals;
handles = guidata(hObject);
viewNum = handles.viewNum;
thisView = viewGet(viewNum,'view');
MLR.clipboard = viewGet(thisView,'baseSform');


% --------------------------------------------------------------------
function pasteBaseSformCallBack(hObject, dump)

mrGlobals;
handles = guidata(hObject);
viewNum = handles.viewNum;
thisView = viewGet(viewNum,'view');
if all(size(MLR.clipboard)==4)
    viewSet(thisView,'baseSform',MLR.clipboard);
else
    mrErrorDlg('(paste base sform) Cannot paste. Clipboard does not contain a valid transformation matrix. Use Edit -> Base Anatomy -> Transforms -> Copy sform.')
end
refreshMLRDisplay(viewNum);


% --------------------------------------------------------------------
function editOverlayCallback(hObject, dump)
handles = guidata(hObject);
viewNum = handles.viewNum;
editOverlayGUI(viewNum);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Analysis Menu Callbacks %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --------------------------------------------------------------------
function applyMotionCompTransformsCallBack(hObject, dump)

handles = guidata(hObject);
viewNum = handles.viewNum;
thisView = viewGet(viewNum,'view');

applyMotionCompTransform(thisView);


%------------------------- glmAnalysisCallback Function ------------------------------%
function glmAnalysisCallback(hObject,dump)
thisView = viewGet(getfield(guidata(hObject),'viewNum'),'view');
glmAnalysis(thisView);


%------------------------- RecomputeAnalysisCallback Function ------------------------------%
% this is to replace the GUI function of old event-related analyses
function recomputeAnalysisCallback(hObject,dump)

viewNum = getfield(guidata(hObject),'viewNum');
thisView = viewGet(viewNum,'view');

if viewGet(thisView,'numAnalyses') > 0
  n = viewGet(thisView,'currentAnalysis');
  groupName = viewGet(thisView,'analysisGroupName',n);
  analysisFunction = viewGet(thisView,'analysisFunction',n);
  guiFunction = viewGet(thisView,'analysisGuiFunction',n);
  switch(guiFunction)
    case {'eventRelatedGUI','eventRelatedGlmGUI'}
      guiFunction = 'glmAnalysisGUI';
      analysisFunction = 'glmAnalysis';
  end
  params = viewGet(thisView,'analysisParams',n);
  % params = guiFunction('groupName',groupName,'dummy',0,'params',params,'thisView',thisView);
  evalstring = ['params = ',guiFunction,'(','''','groupName','''',',groupName,','''','params','''',',params,','''','thisView','''',',thisView);'];
  eval(evalstring);
  % params is empty if GUI cancelled
  if isempty(params)
    disp('(glmAnalysis) Analysis cancelled');
    return
  else
    %check if the group has changed, in which case we need to remove the tseriesFile field so that reconcileParams doesn't get confused
    if isfield(params,'groupName') && ~strcmp(groupName,params.groupName) && isfield(params,'tseriesFile')
       params = rmfield(params,'tseriesFile');
    end
    % thisView = analysisFunction(thisView,params);
    evalstring = ['thisView = ',analysisFunction,'(thisView,params);'];
    eval(evalstring);
    refreshMLRDisplay(viewNum);
  end
else
  mrWarnDlg(sprintf('(mrLoadRetGUI) No analyses loaded'));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Overlay Menu Callbacks %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%------------------------- combineOverlaysCallback Function ------------------------------%
function combineOverlaysCallback(hObject,dump)
thisView = viewGet(getfield(guidata(hObject),'viewNum'),'view');
combineTransformOverlays(thisView);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ROI Menu Callbacks %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --------------------------------------------------------------------
function createSingleVoxelsCallBack(hObject,dump)

viewNum = getfield(guidata(hObject),'viewNum');
thisView = viewGet(viewNum,'view');
[thisView userCancel] = newROI(thisView);
if userCancel,return,end
drawROI(thisView,'single voxels',1);
refreshMLRDisplay(viewNum);


% --------------------------------------------------------------------
function addSingleVoxelsCallBack(hObject,dump)

viewNum = getfield(guidata(hObject),'viewNum');
thisView = viewGet(viewNum,'view');
drawROI(thisView,'single voxels',1);
refreshMLRDisplay(viewNum);


% --------------------------------------------------------------------
function removeSingleVoxelsCallBack(hObject,dump)

viewNum = getfield(guidata(hObject),'viewNum');
thisView = viewGet(viewNum,'view');
drawROI(thisView,'single voxels',0);
refreshMLRDisplay(viewNum);


%------------------------- transformROIsCallback Function ------------------------------%
function transformROIsCallback(hObject,dump)
thisView = viewGet(getfield(guidata(hObject),'viewNum'),'view');
transformROIs(thisView);

