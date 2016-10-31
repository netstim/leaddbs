% mrDispOverlay.m
%
%      usage: [view analysis] = mrDispOverlay(overlay,scanNum,groupNum/analysisStruct,<view>)
%      usage: mrDispOverlay(overlay,scanNum,analysisStructure,<view>)
%         by: justin gardner
%       date: 04/04/07
%    purpose: displays an overlay in MrLoadRet
%
%             For example, to display a map for scan 1, group 2:
%             Where map has dimensions size(map) = [x y s]
% 
%             mrDispOverlay(map,1,2)
%
%             If you want to install maps for several scans, make
%             map into a cell array of maps of the correct
%             dimensions and do
%
%             mrDispOverlay(map,[1 2],2);
%
%             If you want to install multiple maps for a single scan, make
%             a cell array of maps of all the maps you want to associate
%             with the scan and do.
%
%             mrDispOverlay(map,1,2,[],'overlayNames',{'map1name','map2name'});
%           
%             If you want to install the overlay into an existing
%             analysis, rather than a group (say erAnal), pass
%             that analysis instead of the group
% 
%             mrDispOverlay(map,1,erAnal);
%
%             If you want to install the map into an existing
%             view, you can pass that (or empty if you don't want
%             to display at all)
%
%             mrDispOverlay(map,1,erAnal,[]);
%
%             If you just want to save an analysis structure without viewing anything:
%
%             v = newView;
%             map = rand(viewGet(v,'scanDims'));
%             scanNum = 1; groupNum = 1;
%             mrDispOverlay(map,scanNum,groupNum,[],'overlayName=rand','saveName=randAnalysis');
%
%             If you want to return an image array that has the overlay on
%             a specified anatomy, you could do the following
%
%             v = newView;
%             v = loadAnat(v,'jg_left_occipital.hdr');
%             map = rand(viewGet(v,'scanDims'));
%             v = mrDispOverlay(map,1,3,v);
%             img = refreshMLRDisplay(viewGet(v,'viewNum'));
%
%             You can also set optional parameters
%
%             Some parameters you can set
%             'overlayName=name'
%             'cmap',colormap
%             'colormapType=setRangeToMax'
%             'range',[0 1]
%             'clip',[0 1]
%             'interrogator=interrogatorName'
%             'analName=analysisName'
%             'd',d 
%              to set an arbitrary params
%              'params.x',x
%
%             to save instead of display set
%             'saveName=filename'
%
%             e.g.
%             mrDispOverlay(map,1,erAnal,[],'overlayName=myMap');
%
%
function [v, analysis] = mrDispOverlay(overlay,scanNum,groupNum,v,varargin)

% check arguments
if nargin < 3
  help mrDispOverlay
  return
end

% get the arguments
eval(evalargs(varargin));

mrLoadRetViewing = 0;
viewShouldBeDeleted = 0;
% start up a mrLoadRet if we are not passed in a view
% or if we are returning an analysis struct (in this
if ieNotDefined('v')
  if ieNotDefined('saveName') && (nargout < 2)
    v = mrLoadRet;
    mrLoadRetViewing = 1;
  else
    % if we are saving, then don't bring up mrLoadRet
    v = newView;
    viewShouldBeDeleted = 1;
  end
else
  % see if the view is one returned from getMLRView (i.e.
  % it has a figure that needs to be refreshed)
  if ~isempty(viewGet(v,'fignum'))
    mrLoadRetViewing = 1;
  end
end


% if groupNum is actually a structure it means we were passed
% in an analysis strututre
if isstruct(groupNum)
  anal = groupNum;
  groupNum = viewGet(v,'groupNum',anal.groupName);
  v = viewSet(v,'curGroup',groupNum);
  v = viewSet(v,'newAnalysis',anal);
end

% now set to the current view and group
v = viewSet(v,'curGroup',groupNum);
mlrGuiSet(v.viewNum,'scanText',scanNum(1));
groupName = viewGet(v,'groupName',groupNum);

% make the overlay into a cell array, if it is not already.
overlay = cellArray(overlay);

% get the min and max of the overlays
minOverlay = inf;
maxOverlay = -inf;
for i = 1:length(overlay)
  minOverlay = min(minOverlay,min(overlay{i}(:)));
  maxOverlay = max(maxOverlay,max(overlay{i}(:)));
end

% get some default parameters
if ieNotDefined('cmap')
  % colormap is made with a little bit less on the dark end
  cmap = hot(312);
  cmap = cmap(end-255:end,:);
end
if ieNotDefined('colormapType')
  colormapType = 'setRangeToMax';
end
if ieNotDefined('range')
  range = [minOverlay maxOverlay];
end
if ieNotDefined('clip')
  clip = [minOverlay maxOverlay];
end
if ieNotDefined('analName')
  analName = 'mrDispOverlayAnal';
end

% first get an overlayName
if ieNotDefined('overlayName')
  if ieNotDefined('overlayNames')
    overlayName = 'mrDispOverlay';
  else
    overlayName = overlayNames{1};
  end
end
overlayName = fixBadChars(overlayName);

% now look for an interrogator name
if ieNotDefined('interrogator')
  if ieNotDefined('anal')
    interrogator = 'timecoursePlot';
  else
    if isfield(anal,'overlays')
      interrogatorFound = 0;
      for i = 1:length(anal.overlays)
	if strcmp(anal.overlays(i).name,overlayName)
	  interrogator = anal.overlays(i).interrogator;
	end
      end
      if ~interrogatorFound
	interrogator = 'timecoursePlot';
      end
    end
  end
end

% if there is only one scan and multiple overlays,
% it means to install many overlays for the single scan
range = cellArray(range);
clip = cellArray(clip);
cmap = cellArray(cmap);
colormapType = cellArray(colormapType);
if length(scanNum) == 1
  % we have multiple overlays
  numOverlays = length(overlay);
  % get overlay names
  if ieNotDefined('overlayNames') 
    for i = 1:numOverlays
      overlayNames{i} = sprintf('%s%i',overlayName,i);
    end
  end
  % make sure to fill out any missing overlay names
  for i = (length(overlayNames)+1):numOverlays
    overlayNames{i} = sprintf('%s%i',overlayName,i);
  end
  % fill out range
  for i = length(range)+1:numOverlays
    range{i} = range{1};
  end
  % fill out clip
  for i = length(clip)+1:numOverlays
    clip{i} = clip{1};
  end
  % fill out cmap
  for i = length(cmap)+1:numOverlays
    cmap{i} = cmap{1};
  end
  % fill out colormapType
  for i = length(colormapType)+1:numOverlays
    colormapType{i} = colormapType{1};
  end
else
  numOverlays = 1;
end
disp(sprintf('(mrDispOverlay) Num overlays: %i',numOverlays));

for onum = 1:numOverlays
  % create the parameters for the overlay
  dateString = datestr(now);
  % set name, if more than one overlay per scan
  % then number them in order.
  if numOverlays == 1
    o(onum).name = overlayName;
  else
    o(onum).name = overlayNames{onum};
  end
  o(onum).function = '';
  o(onum).groupName = groupName;
  o(onum).reconcileFunction = 'defaultReconcileParams';
  o(onum).data = cell(1,viewGet(v,'nScans'));
  if ieNotDefined('params')
    o(onum).params = [];
  else
    o(onum).params = params;
  end
  o(onum).params.scanNum = [];
  o(onum).date = dateString;
  o(onum).range = range{onum};
  o(onum).clip = clip{onum};
  o(onum).colormap = cmap{onum};
  o(onum).alpha = 1;
  o(onum).colormapType = colormapType{onum};
  o(onum).interrogator = interrogator;
  % now setup the date field
  % if we have multiple overlays, then it means
  % we have a single scan which we are installying
  % multiple overlays (this could be relaxed int he
  % future to have multiple overlays for multiple scans)...
  if numOverlays > 1
    o(onum).data{scanNum(1)} = overlay{onum};
    o(onum).params.scanNum = scanNum(1);
  else
    for i = 1:length(scanNum)
      if length(overlay) >= i
	o(onum).data{scanNum(i)} = overlay{i};
	o(onum).params.scanNum(end+1) = scanNum(i);
	% if we are passed in one overlay but multiple scans
	% then we install the same overlay for all the scans
      elseif length(overlay) == 1
	disp(sprintf('(mrDispOverlay) Only 1 overlay passed in. Installing that overlay for scan %i',scanNum(i)));
	o(onum).data{scanNum(i)} = overlay{1};
	o(onum).params.scanNum(end+1) = scanNum(i);
	% otherwise something is wrong
      else
	disp(sprintf('(mrDispOverlay) Only %i overlays found for %i scans',length(overlay),length(scanNum)));
	return
      end
    end
  end
end

% see if we need to install a bogus analysis
if isempty(viewGet(v,'curAnalysis'))
  a.(analName).name = analName;  % This can be reset by editAnalysisGUI
  a.(analName).type = analName;
  a.(analName).groupName = groupName;
  a.(analName).function = '';
  a.(analName).guiFunction = '';
  a.(analName).params.scanNum = scanNum;
  a.(analName).reconcileFunction = 'defaultReconcileParams';
  % if user didn't pass in an overlay, they are probably
  % using this to install something into the d strucutre
  % instead.
  if ~isempty(overlay)
    a.(analName).overlays = o;
  else
    a.(analName).overlays = [];
  end    
  a.(analName).curOverlay = 1;
  a.(analName).date = dateString;
  if ~ieNotDefined('d')
    a.(analName).d{scanNum} = d;
  end
  v = viewSet(v,'newAnalysis',a.(analName));
% or just install this overlay into the current analysis
    else
  % install overlay
  if ~isempty(overlay)
    for onum = 1:length(o)
      thisOverlay = o(onum);  
      v = viewSet(v,'newOverlay',thisOverlay); %no need to validate, as this is performed by viewSet
    end
  end
  % install d structure
  if ~ieNotDefined('d')
    v = viewSet(v,'newd',d,scanNum);
  end
end

% either save it or refresh the display
if ~ieNotDefined('saveName')
  v = viewSet(v,'analysisName',saveName);
  saveAnalysis(v,saveName);
end

% update mr load ret if it is running
if mrLoadRetViewing
  v = viewSet(v,'curGroup',groupNum);
  v = viewSet(v,'curScan',scanNum(1));
  refreshMLRDisplay(v.viewNum);
else
  % set the group and scan of the view
  if ~isempty(groupNum)
    v = viewSet(v,'curGroup',groupNum(1));
  end
  if ~isempty(scanNum)
    v = viewSet(v,'curScan',scanNum(1));
  end
end

% prepare the output argument
if nargout >= 2
  analysis = viewGet(v,'analysis',viewGet(v,'curAnalysis'));
end

% delete view
if viewShouldBeDeleted
  deleteView(v);
else
end