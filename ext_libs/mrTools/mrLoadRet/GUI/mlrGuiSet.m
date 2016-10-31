function mlrGuiSet(view,field,value,varargin)
%
%        $Id: mlrGuiSet.m 2875 2013-09-27 12:05:35Z julien $
% mlrGuiSet(view,field,value);
%
% view can be either a view structure or a viewNum. Either way, sets
% handles in the global variable: guidata(MLR.views{viewNum}.figure).
% 
% djh 6/2005

mrGlobals;

if ieNotDefined('view'), mrErrorDlg('No view specified.'); end
if ieNotDefined('field'), mrErrorDlg('No parameter specified'); end
if ieNotDefined('value'), val = []; end

% viewNum can be either a view structure or a view number.
if isnumeric(view)
  viewNum = view;
  % *** Delete the following line after eliminating references to view ***
  view = MLR.views{viewNum};
elseif isview(view)
  viewNum = view.viewNum;
else
  return
end
if (viewNum < 1) | (viewNum > length(MLR.views))
  mrErrorDlg('Invalid viewNum,');
end

% Return if no gui
if isempty(MLR.views{viewNum}.figure)
  return
else
  % Get the gui handles to
  handles = guidata(MLR.views{viewNum}.figure);
end

switch lower(field)

 case {'grouppopup','groupstring'}
  % Set the groupPopup string array
  set(handles.groupPopup,'String',value);

 case {'group'}
  % Choose the group
  set(handles.groupPopup,'Value',value);

 case {'roipopup','roistring'}
  % Set the roiPopup string array
  set(handles.roiPopup,'String',value);

 case {'roi'}
  % Choose the roi
  if strcmp(get(handles.roiPopup,'style'),'popupmenu')
    value = value(1); %if this is not a listbox, we set only one ROI;
  end
  set(handles.roiPopup,'Value',value);

 case {'basepopup'}
  % Set the basePopup string array
  set(handles.basePopup,'String',value);
  % see if there are any registered callbacks
  if ~isempty(viewGet(view,'figNum'))
    callbacks = viewGet(view,'callback','baseChange');
    % and call them
    for iCallback = 1:length(callbacks)
      view = feval(callbacks{iCallback},view);
    end
  end
 case {'labelrois'}
  % mlrGuiSet(view,'labelrois',value);
  if value
    set(handles.labelsROIsMenuItem,'Checked','on');
  else
    set(handles.labelsROIsMenuItem,'Checked','off');
  end
  if isfield(handles,'displayROILabels')
    set(handles.displayROILabels,'value',value);
  end
  
 case {'showrois'}
  % mlrGuiSet(view,'showrois',value);
  
  % figure out which menu item should be checked
  onItem = find(strcmp(value,{'all','all perimeter','selected','selected perimeter','group','group perimeter','hide'}));
  onOrOff = {'off','off','off','off','off','off','off'};
  onOrOff{onItem} = 'on';
  % turn the check marks on/off
  set(handles.showAllMenuItem,'Checked',onOrOff{1});
  set(handles.showAllPerimeterMenuItem,'Checked',onOrOff{2});
  set(handles.showSelectedMenuItem,'Checked',onOrOff{3});
  set(handles.showSelectedPerimeterMenuItem,'Checked',onOrOff{4});
  set(handles.showGroupMenuItem,'Checked',onOrOff{5});
  set(handles.showGroupPerimeterMenuItem,'Checked',onOrOff{6});
  set(handles.hideROIsMenuItem,'Checked',onOrOff{7});
  
  if isfield(handles,'roiDisplayModePopup')
    switch(value)
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
    set(handles.roiDisplayModePopup,'value',roiDisplayMode);
  end
  
 case {'basetype'}
  % mlrGuiSet(view,'baseType',value);
  % value = 0 for regular or 1 for flat, 2 for surface
  if value == 0
    set(handles.sliceText,'Visible','on');
    set(handles.sliceSlider,'Visible','on');
    set(handles.slice,'Visible','on');    
    set(handles.sagittalRadioButton,'Visible','on');
    set(handles.coronalRadioButton,'Visible','on');
    set(handles.axialRadioButton,'Visible','on');
    set(handles.corticalDepth,'Visible','off');
    set(handles.corticalDepthSlider,'Visible','off');
    set(handles.corticalDepthText,'Visible','off');
    set(handles.flatViewerMenuItem,'Enable','off');
    set(handles.calcDistMenu, 'Enable', 'off');
    set(handles.convertCorticalDepthRoiMenuItem,'Enable','off');
    set(handles.axisSingle,'Visible','on');
    set(handles.axisMulti,'Visible','on');
    set(handles.axis3D,'Visible','on');
    if isfield(handles,'corticalMaxDepthSlider')
      set(handles.corticalMaxDepthSlider,'Visible','off');
      set(handles.corticalMaxDepthText,'Visible','off');
      set(handles.linkMinMaxDepthCheck,'Visible','off');
    end
  elseif value >= 1
    set(handles.sliceText,'Visible','off');
    set(handles.sliceSlider,'Visible','off');
    set(handles.slice,'Visible','off');    
    set(handles.sagittalRadioButton,'Visible','off');
    set(handles.coronalRadioButton,'Visible','off');
    set(handles.axialRadioButton,'Visible','off');
    set(handles.corticalDepth,'Visible','on');
    set(handles.corticalDepthSlider,'Visible','on');
    set(handles.corticalDepthText,'Visible','on');
    set(handles.flatViewerMenuItem,'Enable','on');
    set(handles.flatViewerMenuItem,'Label','Flat Viewer');
    set(handles.calcDistMenu, 'Enable', 'on');
    set(handles.convertCorticalDepthRoiMenuItem,'Enable','on');
    set(handles.axisSingle,'Visible','off');
    set(handles.axisMulti,'Visible','off');
    set(handles.axis3D,'Visible','off');
    if isfield(handles,'corticalMaxDepthSlider')
      set(handles.corticalMaxDepthSlider,'Visible','on');
      set(handles.corticalMaxDepthText,'Visible','on');
      set(handles.linkMinMaxDepthCheck,'Visible','on');
    end
  end
  if value == 2
    set(handles.createLineMenuItem,'Enable','off');
    set(handles.createRectangleMenuItem,'Enable','off');
    set(handles.createContiguousMenuItem,'Enable','off');
    set(handles.addLineMenuItem,'Enable','off');
    set(handles.addRectangleMenuItem,'Enable','off');
    set(handles.addContiguousMenuItem,'Enable','off');
    set(handles.removeLineMenuItem,'Enable','off');
    set(handles.removeRectangleMenuItem,'Enable','off');
    set(handles.removeContiguousMenuItem,'Enable','off');
    if isfield(handles,'createSingleVoxelsRoiMenuItem')
      set(handles.createSingleVoxelsRoiMenuItem,'Enable','off')
      set(handles.addSingleVoxelsRoiMenuItem,'Enable','off')
      set(handles.removeSingleVoxelsRoiMenuItem,'Enable','off')
    end
    set(handles.rotateSlider,'SliderStep',[15 45]./360);
    set(handles.baseTiltSlider,'SliderStep',[15 45]./360);
    set(handles.baseTiltSlider,'Visible','on');
    set(handles.baseTiltText,'Visible','on');
    set(handles.baseTilt,'Visible','on');
    set(handles.flatViewerMenuItem,'Enable','on');
    set(handles.flatViewerMenuItem,'Label','Surface Viewer');
    set(handles.calcDistMenu, 'Enable', 'off');
    set(handles.convertCorticalDepthRoiMenuItem,'Enable','on');
  else
    set(handles.baseTiltSlider,'Visible','off');
    set(handles.baseTiltText,'Visible','off');
    set(handles.baseTilt,'Visible','off');
    set(handles.createLineMenuItem,'Enable','on');
    set(handles.createRectangleMenuItem,'Enable','on');
    set(handles.createContiguousMenuItem,'Enable','on');
    set(handles.addLineMenuItem,'Enable','on');
    set(handles.addRectangleMenuItem,'Enable','on');
    set(handles.addContiguousMenuItem,'Enable','on');
    set(handles.removeLineMenuItem,'Enable','on');
    set(handles.removeRectangleMenuItem,'Enable','on');
    set(handles.removeContiguousMenuItem,'Enable','on');
    if isfield(handles,'createSingleVoxelsRoiMenuItem')
      set(handles.createSingleVoxelsRoiMenuItem,'Enable','on')
      set(handles.addSingleVoxelsRoiMenuItem,'Enable','on')
      set(handles.removeSingleVoxelsRoiMenuItem,'Enable','on')
    end
    set(handles.rotateSlider,'SliderStep',[1 45]./360);
  end		  
 case {'basevolume'}
  % Choose the baseVolume
  set(handles.basePopup,'Value',value);
  
 case {'basedims'}
  % mlrGuiSet(view,'baseDims',[ydim xdim zdim]);
  newDims = value;
  newCoords = min(handles.coords,newDims);
  handles.coords = min(handles.coords,newCoords);

 case {'curcoords'}
  % mlrGuiSet(view,'curcoords',[ydim xdim zdim]);
  if length(value)==3
    handles.coords = value;
  end

 case {'basegamma'}
  % mlrGuiSet(view,'baseGamma',value);
  set(handles.baseGammaSlider,'Value',value);
  set(handles.baseGammaText,'String',thisNum2str(value));
  
 case {'basetilt'}
  % mlrGuiSet(view,'baseMax',value);
  set(handles.baseTiltSlider,'Value',value);
  set(handles.baseTiltText,'String',thisNum2str(value));

 case {'analysispopup'}
  % mlrGuiSet(view,'analysisPopup',strings);
  set(handles.analysisPopup,'String',value);

 case {'analysis'}
  % mlrGuiSet(view,'analysis',overlayNum);
  set(handles.analysisPopup,'Value',value);
  
 case {'overlaypopup'}
  % mlrGuiSet(view,'overlayPopup',strings);
  % mlrGuiSet(view,'overlayPopup',strings,overlayList); 
  %optional argument overlayList indicates that this is a subset of strings to change
  if ~strcmp(value,'none') 
    if ieNotDefined('varargin')
% %       overlayList = 1:length(value);
    else
      overlayList = varargin{1};
      newStrings=  value;
      value = get(handles.overlayPopup,'String');
      value(overlayList) = newStrings;
    end
% %     %identify overlays that have been masked by putting a star before their name
% %     epsilon = 1e-7; %value differing by less than epsilon are considered equal
% %     for iOverlay = overlayList
% %       clip = viewGet(view,'overlayclip',iOverlay);
% %       minOverlayData = viewGet(view,'minoverlaydata',iOverlay);
% %       maxOverlayData = viewGet(view,'maxoverlaydata',iOverlay);
% %       if (~isempty(minOverlayData) && (clip(1)-minOverlayData)>epsilon) ||...
% %             (~isempty(maxOverlayData) && (maxOverlayData-clip(2))>epsilon) || ...
% %             clip(1)==clip(2) %if min and max clip values are equal, the whole overlay will be masked
% %          value{iOverlay} = [char(42) ' ' value{iOverlay}];
% %       else
% %          value{iOverlay} = ['  ' value{iOverlay}];
% %       end
% %     end
  else
    set(handles.overlayPopup,'value',1);
  end
  % for matlab version 2014a and above, the listboxtop property is not
  % correctly updated when setting the value or the strings
  if ~verLessThan('matlab','8.3') && strcmp(get(handles.overlayPopup,'style'),'listbox')
    set(handles.overlayPopup,'ListboxTop',1) %need to set it manually otherwise a warning will be issued
  end
  set(handles.overlayPopup,'String',value);

 case {'overlay'}
  % mlrGuiSet(view,'overlay',overlayNum);
  if strcmp(get(handles.overlayPopup,'style'),'popupmenu')
    value = value(1); %if this is not a listbox, we set only one overlay;
  end
  set(handles.overlayPopup,'Value',value);
    
  if length(value)==1
    set(handles.overlayMinSlider,'enable','on')
    set(handles.overlayMinText,'enable','on')
    set(handles.overlayMaxSlider,'enable','on')
    set(handles.overlayMaxText,'enable','on')
    if isfield(handles,'clippingOverlaysListbox')
      set(handles.clippingOverlaysListbox,'enable','on');
    end
%     set(handles.alphaText,'enable','on');
%     set(handles.alphaSlider,'enable','on');
  else %if more than one overlay selected, disable overlay controls
    if isfield(handles,'clippingOverlaysListbox')
      set(handles.clippingOverlaysListbox,'enable','off');
%       set(handles.clippingOverlaysListbox,'String',[],'enable','off');
    end
    set(handles.overlayMinSlider,'enable','off')
    set(handles.overlayMinText,'enable','off')
    set(handles.overlayMaxSlider,'enable','off')
    set(handles.overlayMaxText,'enable','off')
%     set(handles.alphaText,'enable','off');
%     set(handles.alphaSlider,'enable','off');
  end
  
 case {'overlaymin'}
  % mlrGuiSet(view,'overlayMin',value);
  if rem(value,1)~=0 %if the value is not an integer
    value = floor(double(value)*1e6)/1e6; %round it down
  end
  value = clipToSlider(handles.overlayMinSlider,value);
  set(handles.overlayMinSlider,'Value',value);
  set(handles.overlayMinText,'String',thisNum2str(value)); 

 case {'overlayminrange'}
  % mlrGuiSet(view,'overlayMinRange',[min,max]);
  if ~isempty(value)
    if ~all(isfinite(value))
      mrWarnDlg(sprintf('(mlrGuiSet) Cannot display Overlay Min Slider because overlay range is not finite [%f %f]',value(1),value(2)));
      set(handles.overlayMinSlider,'Min',max(value(1),-realmax(class(value))),'Max',min(value(2),realmax(class(value))),'visible','off');
    elseif value(2) < value(1)
      mrWarnDlg(sprintf('(mlrGuiSet) Cannot display Overlay Min Slider because overlay range is not increasing [%f > %f]',value(1),value(2)));
      set(handles.overlayMinSlider,'Min',max(value(1),-realmax(class(value))),'Max',min(value(2),realmax(class(value))),'visible','off');
    else
      set(handles.overlayMinSlider,'Min',value(1),'Max',value(2),'visible','on');
      end
  end
 case {'overlaymax'}
  % mlrGuiSet(view,'overlayMax',value);
  if ~isempty(value)
    if rem(value,1)~=0 %if the value is not an integer
      value = ceil(double(value)*1e6)/1e6; %round it up
    end
    value = clipToSlider(handles.overlayMaxSlider,value);
    set(handles.overlayMaxSlider,'Value',value);
    set(handles.overlayMaxText,'String',thisNum2str(value)); 
  end
 case {'overlaymaxrange'}
  % mlrGuiSet(view,'overlayMinRange',[min,max]);
  if ~isempty(value)
    if ~all(isfinite(value))
      mrWarnDlg(sprintf('(mlrGuiSet) Cannot display Overlay Max Slider because overlay range is not finite [%f %f]',value(1),value(2)));
      set(handles.overlayMaxSlider,'Min',max(value(1),-realmax(class(value))),'Max',min(value(2),realmax(class(value))),'visible','off');
    elseif value(2) < value(1)
      mrWarnDlg(sprintf('(mlrGuiSet) Cannot display Overlay Max Slider because overlay range is not increasing [%f > %f]',value(1),value(2)));
      set(handles.overlayMaxSlider,'Min',max(value(1),-realmax(class(value))),'Max',min(value(2),realmax(class(value))),'visible','off');
    else
      set(handles.overlayMaxSlider,'Min',value(1),'Max',value(2),'visible','on');
    end
  end
 case {'clipacrossoverlays'}
  % mlrGuiSet(view,'clipAcrossOverlays',value);
  if isfield(handles,'clipAcrossOverlays') 
    set(handles.clipAcrossOverlays,'value',value)
  end

 case {'multiaxis'}
  % mlrGuiSet(view,'multiAxis',value);

  % set the radio buttons appropriately
  controlNames = {'Single','Multi','3D'};
  for iControl = 1:3
    if iControl~=(value+1)
      set(handles.(['axis' controlNames{iControl}]),'Value',0);
    else
      set(handles.(['axis' controlNames{iControl}]),'Value',1);
    end
  end

  
  if value == 1
    % create the locations for the 3 axis
    handles.anatMultiPosition(3,:) = handles.anatPosition;
    handles.anatMultiPosition(3,3) = (handles.anatMultiPosition(3,3)-handles.marginSize)/2;
    handles.anatMultiPosition(3,4) = (handles.anatMultiPosition(3,4)-handles.marginSize)/2;
    handles.anatMultiPosition(2,:) = handles.anatMultiPosition(3,:);
    handles.anatMultiPosition(2,2) = handles.anatMultiPosition(2,2)+(handles.anatPosition(4)+handles.marginSize)/2;
    handles.anatMultiPosition(1,:) = handles.anatMultiPosition(2,:);
    handles.anatMultiPosition(1,1) = handles.anatMultiPosition(1,1)+(handles.anatPosition(3)+handles.marginSize)/2;
    handles.anatMultiPosition(4,:) = handles.anatMultiPosition(3,:);
    handles.anatMultiPosition(4,1) = handles.anatMultiPosition(1,1);
    % set the regular window to the bottom right
    set(handles.axis,'Position',handles.anatMultiPosition(4,:));
    % display the axis and clear
    set(handles.sliceAxis(1),'Position',handles.anatMultiPosition(1,:));set(handles.sliceAxis(1),'Visible','on');cla(handles.sliceAxis(1));axis(handles.sliceAxis(1),'off');
    set(handles.sliceAxis(2),'Position',handles.anatMultiPosition(2,:));set(handles.sliceAxis(2),'Visible','on');cla(handles.sliceAxis(2));axis(handles.sliceAxis(2),'off');
    set(handles.sliceAxis(3),'Position',handles.anatMultiPosition(3,:));set(handles.sliceAxis(3),'Visible','on');cla(handles.sliceAxis(3));axis(handles.sliceAxis(3),'off');
    % turn on tilt slider
    set(handles.rotateSlider,'SliderStep',[15 45]./360);
    set(handles.baseTiltSlider,'SliderStep',[15 45]./360);
    set(handles.baseTiltSlider,'Visible','on');
    set(handles.baseTiltText,'Visible','on');
    set(handles.baseTilt,'Visible','on');
  else
    % set the regular window to the normal position
    set(handles.axis,'Position',handles.anatPosition);
    cla(handles.axis,'reset');
    axis(handles.axis,'off');
    % clear and hide the other axes
    cla(handles.sliceAxis(1),'reset');set(handles.sliceAxis(1),'Visible','off');
    cla(handles.sliceAxis(2),'reset');set(handles.sliceAxis(2),'Visible','off');
    cla(handles.sliceAxis(3),'reset');set(handles.sliceAxis(3),'Visible','off');
    if value == 2
      % turn on tilt slider
      set(handles.rotateSlider,'SliderStep',[15 45]./360);
      set(handles.baseTiltSlider,'SliderStep',[15 45]./360);
      set(handles.baseTiltSlider,'Visible','on');
      set(handles.baseTiltText,'Visible','on');
      set(handles.baseTilt,'Visible','on');
    else
      % turn off tilt slider
      set(handles.rotateSlider,'SliderStep',[1 45]./360);
      %set(handles.baseTiltSlider,'Visible','off');
      %set(handles.baseTiltText,'Visible','off');
      %set(handles.baseTilt,'Visible','off');
    end
  end
 case {'clippingoverlays'}
  % mlrGuiSet(view,'clippingOverlays',overlayList);
  overlayList=value;
  if isfield(handles,'clippingOverlaysListbox') 
    overlayStrings=get(handles.overlayPopup,'String');
    selected=get(handles.clippingOverlaysListbox,'value');
    if selected>length(overlayList)
      selected=1;
    end
    if viewGet(view,'nOverlays')
      %identify overlays that have been masked by putting a star before their name
      epsilon = 1e-7; %value differing by less than epsilon are considered equal
      for iOverlay = overlayList
        clip = viewGet(view,'overlayclip',iOverlay);
        minOverlayData = viewGet(view,'minoverlaydata',iOverlay);
        maxOverlayData = viewGet(view,'maxoverlaydata',iOverlay);
        if (~isempty(minOverlayData) && (clip(1)-minOverlayData)>epsilon) ||...
              (~isempty(maxOverlayData) && (maxOverlayData-clip(2))>epsilon) || ...
              (~isempty(clip) && clip(1)==clip(2)) %if min and max clip values are equal, the whole overlay will be masked
           overlayStrings{iOverlay} = [char(42) ' ' overlayStrings{iOverlay}];
        else
           overlayStrings{iOverlay} = ['  ' overlayStrings{iOverlay}];
        end
      end
    end
    set(handles.clippingOverlaysListbox,'String',overlayStrings(value),'value',selected);
  end

 case {'alpha','overlayalpha'}
  % mlrGuiSet(view,'alpha',value);
  value = clipToSlider(handles.alphaSlider,value);
  set(handles.alphaSlider,'Value',value);
  set(handles.alphaText,'String',thisNum2str(value));
  set(handles.alphaSlider,'sliderStep',[0.1 0.5]);

 case {'addpanel'}
    % add a panel - this is to display a set of controls on the
    % right - this is used by mlrAdjustGUI
    % mlrGuiSet(v,'addPanel','panelName',.5)
    % where panelName is an arbitraty name of panel
    % .5 is the percent height of the panel (has to be between 0 and 1)
    if length(varargin)<1
      disp(sprint('(mlrGuiSet) addPanel must specify percent height'));
      return
    end
    % check argument for panelSize
    if length(varargin)>=1
      percentSize = varargin{1};
      if ~isnumeric(percentSize) || (length(percentSize)~=1) || (percentSize<0) || (percentSize>1)
	disp(sprintf('(mlrGuiSet) Panel must specify a percent size between 0 and 1'));
	return
      end
    else
      percentSize = 0.25;
    end
    % create a panel item
    f = viewGet(view,'figNum');
    if ~isempty(f)
      panelHandle = uipanel(f,'Visible','off');
      % add the fields to the global
      if ~isfield(MLR,'panels')
	MLR.panels{1} = {value,panelHandle,percentSize,false};
      else
	MLR.panels{end+1} = {value,panelHandle,percentSize,false};
      end
      % display the panel
      mlrGuiSet(view,'dispPanel',value);
      handles = guidata(viewGet(view,'figNum'));
    end
      
 case {'disppanel'}
    % display a panel - this is used to display a panel that has
    % been added by addPanel
    panelY = 1;
    % make the axis display shorter
    handles.anatPosition(3) = 1-handles.anatPosition(1)-0.25-handles.marginSize;
    % save the anat position
    guidata(viewGet(view,'figNum'),handles);
    % make the colorbar display shorter
    colorbarPosition = get(handles.colorbar,'Position');
    colorbarPosition(3) = 1-colorbarPosition(1)-0.25-handles.marginSize;
    set(handles.colorbar,'Position',colorbarPosition);
    if isfield(handles,'colorbarRightBorder')
      widthBorder = 0.001;
      borderPosition = colorbarPosition;
      borderPosition(1) = borderPosition(1)+borderPosition(3)-widthBorder;
      borderPosition(3) = widthBorder;
      set(handles.colorbarRightBorder,'Position',borderPosition);
    end
    refreshMLRDisplay(viewGet(view,'viewNum'));
    for iPanel = 1:length(MLR.panels)
      % if panel is visible, then shift panelY accordingly
      if ~strcmp(MLR.panels{iPanel}{1},value) && MLR.panels{iPanel}{4}
	panelY = panelY - MLR.panels{iPanel}{3};
	% found panel to display
      elseif strcmp(MLR.panels{iPanel}{1},value) 
	if (panelY - MLR.panels{iPanel}{3}) >= 0
	  set(MLR.panels{iPanel}{2},'Position',[.75 (panelY-MLR.panels{iPanel}{3}-handles.marginSize) .25-handles.marginSize (1-2*handles.marginSize)*MLR.panels{iPanel}{3}]);
	  set(MLR.panels{iPanel}{2},'Visible','on');
	  % set the field so that we know this one is visible
	  MLR.panels{iPanel}{4} = true;
	  % FIX, FIX, FIX - need to shift down panels below this one
	else
	  disp(sprintf('(mlrGuiSet) Could not display panel %s, because there is not enough space',value));
	end
      end
    end
 case {'hidepanel'}
    % hide panel - this is used to hide a panel that has
    % been added by addPanel
    panelsDisplaying = 0;panelFound = false;
    if isfield(MLR,'panels')
      for iPanel = 1:length(MLR.panels)
	% make panel invisible if we found a match
	if strcmp(MLR.panels{iPanel}{1},value) 
	  set(MLR.panels{iPanel}{2},'Visible','off');
	  % set the field to say that it is not displaying
	  MLR.panels{iPanel}{4} = false;
	  % we have found our panel
	  panelFound = true;
	  % not the panel to hide and is being displayed
	elseif MLR.panels{iPanel}{4}
	  panelsDisplaying = true;
	  % FIX, FIX, FIX - need to move up panels below this one
	  if panelFound
	  end
	end
      end
    end
    % if there are no panels displaying, then
    % can show the anatomy and colorbar wider
    % make the axis display shorter
    handles.anatPosition(3) = 1-handles.anatPosition(1)-handles.marginSize;
    % save back anatPosition
    guidata(viewGet(view,'figNum'),handles);
    % make the colorbar display shorter
    colorbarPosition = get(handles.colorbar,'Position');
    if ~isfield(handles,'colorbarRightBorder')
      colorbarPosition(3) = 1-colorbarPosition(1)-handles.marginSize;
    else
      %if there is a right colorscale border, this means that multiple colorbar will be displayed 
      %and we need a bit more space on the right for the scale value 
      colorbarPosition(3) = 1-colorbarPosition(1)-0.05-handles.marginSize;
      %now change the position of the right colorscale border
      widthBorder = 0.001;
      borderPosition = colorbarPosition;
      borderPosition(1) = borderPosition(1)+borderPosition(3)-widthBorder;
      borderPosition(3) = widthBorder;
      set(handles.colorbarRightBorder,'Position',borderPosition);
    end
    set(handles.colorbar,'Position',colorbarPosition);
    refreshMLRDisplay(viewGet(view,'viewNum'));
 case {'nscans'}
  % mlrGuiSet(view,'nscans',value);
  nScans = round(value);
  curScan = round(str2num(get(handles.scanText,'String')));
  if (nScans > 1)
    set(handles.scanSlider,'Min',1);
    set(handles.scanSlider,'Max',nScans);
    set(handles.scanSlider,'sliderStep',[1/(nScans-1) 1/(nScans-1)]);          
    set(handles.scanSlider,'Visible','on');
    curScan = min(curScan,nScans);
  else
    set(handles.scanSlider,'Min',0.9);
    set(handles.scanSlider,'Max',1.1);
    set(handles.scanSlider,'Value',1); %this wasn't needed until matlab ver 2014a, but now, visible[off] doesn't seem to work when 'value' is not within 'range' and therefore the slider is not rendered 
    set(handles.scanSlider,'Visible','off');
  end
  mlrGuiSet(view,'scan',curScan);

 case {'scan'}
  % mlrGuiSet(view,'scan',value);
  set(handles.scanSlider,'Value',clipToSlider(handles.scanSlider,value,1));
  set(handles.scanText,'String',num2str(value));
  % description
  description = viewGet(view,'description',value);
  set(viewGet(view,'figNum'),'Name',sprintf('%s: %s',getLastDir(MLR.homeDir),description));

 case {'scantext'}
  % mlrGuiSet(view,'scanText',value);
  set(handles.scanText,'String',num2str(value));
  % description
  description = viewGet(view,'description',value);
  set(viewGet(view,'figNum'),'Name',sprintf('%s: %s',getLastDir(MLR.homeDir),description));
  
 case {'nslices'}
  % mlrGuiSet(view,'nslices',value);
  value = round(value);
  if (value > 1)
    % reset range
    set(handles.sliceSlider,'Min',1);
    set(handles.sliceSlider,'Max',value);
    % reset step
    set(handles.sliceSlider,'sliderStep',[1/(value-1),1/(value-1)]);
    set(handles.sliceSlider,'Visible','on');
  else
    set(handles.sliceSlider,'Visible','off');
    set(handles.sliceText,'String',0);
  end

 case {'slice'}
  % mlrGuiSet(view,'slice',value);
  value = clipToSlider(handles.sliceSlider,value,1);
  if ~isempty(value)
    handles.coords(handles.sliceOrientation) = value;
    set(handles.sliceSlider,'Value',value);
    set(handles.sliceText,'String',num2str(value));
  end
 case {'slicetext'}
  % mlrGuiSet(view,'sliceText',value);
  handles.coords(handles.sliceOrientation) = value;
  set(handles.sliceText,'String',num2str(value));
  view = viewSet(view,'curSlice',value);

 case {'corticaldepth'} %this sets a single cortical depth, whether there are one or two sliders
  % mlrGuiSet(view,'corticalDepth',value);
   mlrGuiSet(view,'corticalMinDepth',value);
   if isfield(handles,'linkMinMaxDepthCheck') 
     set(handles.linkMinMaxDepthCheck,'value',true)
   end
   if isfield(handles,'corticalMaxDepthSlider')
     mlrGuiSet(view,'corticalMaxDepth',value);
   end
  
 case {'corticalmindepth'}
  % mlrGuiSet(view,'corticalDepth',value);
  value = min(value,1);value = max(value,0);
  set(handles.corticalDepthSlider,'Value',value);
  set(handles.corticalDepthText,'String',num2str(value));
  
 case {'corticalmaxdepth'}
  % mlrGuiSet(view,'corticalDepth',value);
  if isfield(handles,'corticalMaxDepthSlider')
    value = min(value,1);value = max(value,0);
    set(handles.corticalMaxDepthSlider,'Value',value);
    set(handles.corticalMaxDepthText,'String',num2str(value));
  end
  
 case {'sliceorientation'}
  % mlrGuiSet(view,'sliceorientation',value);
  sliceOrientation = value;
  handles.sliceOrientation = sliceOrientation;
  % Set the correct radio button and unset the other two
  switch sliceOrientation
    % axial
   case 1
    % sagittal
    set(handles.sagittalRadioButton,'Value',1);
    set(handles.coronalRadioButton,'Value',0);
    set(handles.axialRadioButton,'Value',0);
   case 2
    % coronal
    set(handles.sagittalRadioButton,'Value',0);
    set(handles.coronalRadioButton,'Value',1);
    set(handles.axialRadioButton,'Value',0);
   case 3
    % axial
    set(handles.sagittalRadioButton,'Value',0);
    set(handles.coronalRadioButton,'Value',0);
    set(handles.axialRadioButton,'Value',1);
  end
  
 case {'rotate'}
  % mlrGuiSet(view,'rotate',value);
  
  value = clipToSlider(handles.rotateSlider,value);
  set(handles.rotateText,'String',num2str(value));
  set(handles.rotateSlider,'Value',value);

 case {'viewnum'}
  % mlrGuiSet(view,'viewnum',value);
  handles.viewNum = value;

 otherwise
  error(['Invalid field: ',field]);
end
guidata(MLR.views{viewNum}.figure,handles);


function value = clipToSlider(slider,value,integerFlag)
% Clips value so that it doesn' texceed slider limits.
% slider is a slider handle
% value must be a number (otherwise, use current slider value)
% integerFlag forces value to be an integer
if ieNotDefined('integerFlag')
  integerFlag = 0;
end
if ~isnumeric(value)
  value = get(slider,'Value');
end
if integerFlag
  if (value < get(slider,'Min'))
    value = ceil(get(slider,'Min'));
  end
  if (value > get(slider,'Max'))
    value = floor(get(slider,'Max'));
  end
else
  if (value < get(slider,'Min'))
    value = get(slider,'Min');
  end
  if (value > get(slider,'Max'))
    value = get(slider,'Max');
  end
end

%modified num2str to increase the number of decimals for reals
function value = thisNum2str(value)

  if rem(value,1)~=0
    value = num2str(value,'%.6f');
  else
    value = num2str(value);
  end

