% searchlightClassification.m
%
%        $Id: searchlightClassification.m 1839 2010-11-14 17:45:36Z julien $
%      usage: thisView = eventRelated(thisView,params)
%         by: alex beckett
%       date: 10/20/06
%    purpose: event related data analysis
%
%             if you just want a default parameter structure you
%             can do:
% 
%             v = newView;
%             [v params] = searchlightClassification(v,[],'justGetParams=1','defaultParams=1','scanList=1')
%
%             Note that justGetParams,defualtParams and scanList are independent parameters, so
%             if you want, say to bring up the GUI to set the params, but not run the analysis, you
%             can do:
%             [v params] = searchlightClassification(v,[],'justGetParams=1');
%
function [thisView d] = corticalSearchlight(thisView,params,varargin)
% check arguments
if ~any(nargin == [1 2 3 4 5])
  help searchlightClassification
  return
end

mrGlobals;

% other arguments
eval(evalargs(varargin));%,[],[],{'justGetParams','defaultParams','scanList'}));
if ieNotDefined('justGetParams'),justGetParams = 0;end
if ieNotDefined('defaultParams'),defaultParams = 0;end
if ieNotDefined('scanList'),scanList = [];end
if ieNotDefined('params'),params = [];end

% First get parameters
if isempty(params) || justGetParams
    params = searchlightClassGUI('thisView',thisView,'params',params,'defaultParams',defaultParams,'scanList',scanList);
%   else
%     params = searchlightClassGUI('groupName',viewGet(thisView,'groupName'),'scanList',scanList,'roilist',viewGet(thisView,'roiNames'));
%   end
end

% Abort if params empty
if ieNotDefined('params')
  disp('(searchlightAnalysis) Searchlight Analysis cancelled');
  return
% just return parameters
elseif justGetParams
%   return
end


% set the group
thisView = viewSet(thisView,'groupName',params.groupName);
% Reconcile params with current status of group and ensure that it has
% the required fields. 
params = defaultReconcileParams([],params);

params.fweMethod='Adaptive Step-down';
params.fdrAssumption= 'Independence/Positive dependence';
params.fdrMethod= 'Adaptive Step-up';
% params.fweAdjustment= 1;
% params.fdrAdjustment =  1;
params.testOutput= 'Z value';
params.trueNullsEstimationMethod= 'Least Squares';
params.trueNullsEstimationThreshold= 0.0500;

set(viewGet(thisView,'figNum'),'Pointer','watch');drawnow;
for scanNum = params.scanNum
  d = loadScan(thisView,scanNum,[],0);
  d = getStimvol(d,params.scanParams{scanNum});
  % do any called for preprocessing
  d = eventRelatedPreProcess(d,params.scanParams{scanNum}.preprocess);

  [out{scanNum},params]=classify_corticalSearchlight(thisView,d,params,params.scanParams{scanNum});
  if isempty(out{scanNum})
    disp('Cortical Searchlight Analysis cancelled');
    set(viewGet(thisView,'figNum'),'Pointer','arrow');
    return
  end
end

%-------------------------------------------------------- Output Analysis ---------------------------------------------------
dateString = datestr(now);
classAnal.name = params.saveName;
classAnal.type = 'searchClass';
classAnal.groupName = params.groupName;
classAnal.function = 'corticalSearchlight';
classAnal.reconcileFunction = 'defaultReconcileParams';
classAnal.mergeFunction = 'defaultMergeParams';
classAnal.guiFunction = 'searchlightClassGUI';
classAnal.params = params;
classAnal.date = dateString;


%--------------------------------------------------------- Output overlay structures
nScans = viewGet(thisView,'nScans');
% create generic parameters 
defaultOverlay.groupName = params.groupName;
defaultOverlay.function = 'corticalSearchlight';
defaultOverlay.reconcileFunction = 'defaultReconcileParams';
defaultOverlay.date = dateString;
defaultOverlay.params = cell(1,nScans);
% colormap is made with a little bit less on the dark end
defaultOverlay.colormap = hot(312);
defaultOverlay.colormap = defaultOverlay.colormap(end-255:end,:);
defaultOverlay.alpha = 1;
defaultOverlay.interrogator = 'timecoursePlot';
defaultOverlay.mergeFunction = 'defaultMergeParams';
defaultOverlay.colormapType = 'normal';
defaultOverlay.range = [0 1];
defaultOverlay.clip = [0 1];
defaultOverlay.alphaOverlay='';
defaultOverlay.alphaOverlayExponent=1;
defaultOverlay.data = cell(1,nScans);
defaultOverlay.name = '';
% for iScan = params.scanNum
%    defaultOverlay.data{iScan} = NaN(scanDims); %to make values outside the box transparent
% end

%------------------------------------save overlays
overlays = defaultOverlay;
overlays.name = ['mean accuracy (radius = ',num2str(params.radius),')'];
for iScan = params.scanNum
   overlays.data{iScan}=out{iScan}.accDiag;
   overlays.params{iScan} = params.scanParams{iScan};
end

thisOverlay = defaultOverlay;
for i=1:params.numberEvents
  overlays(end+1)=thisOverlay;
  overlays(end).name=[params.EVnames{i},'_acc (radius = ',num2str(params.radius),')'];
  for iScan = params.scanNum
      overlays(end).data{iScan}=out{iScan}.accDiag_Class(:,:,:,i);
      overlays(end).params=params.scanParams{iScan};    %JB: not too sure about this 
  end
end

if params.sigTest
    thisOverlay=defaultOverlay;
%     thisOverlay.colormap = statsColorMap(256);
    thisOverlay.range = [0 9];
    thisOverlay.clip = [0 9];
    overlays(end+1)=thisOverlay;
    overlays(end).name=['Z (radius = ',num2str(params.radius),')'];
    for iScan = params.scanNum
        overlays(end).data{iScan}=out{iScan}.pDiag;
        overlays(end).params=params.scanParams{iScan};
    end
    for i=1:params.numberEvents
        overlays(end+1)=thisOverlay;
        overlays(end).name=[params.EVnames{i},'_Z (radius = ',num2str(params.radius),')'];
        for iScan = params.scanNum
            overlays(end).data{iScan}=out{iScan}.pDiag_Class(:,:,:,i);
            overlays(end).params=params.scanParams{iScan};
        end
    end
    if params.fdrAdjustment
        overlays(end+1)=thisOverlay;
        overlays(end).name=['fdr Z (radius = ',num2str(params.radius),')'];
        for iScan = params.scanNum
            overlays(end).data{iScan}=out{iScan}.fdrDiag;
            overlays(end).params=params.scanParams{iScan};
        end
        for i=1:params.numberEvents
            overlays(end+1)=thisOverlay;
            overlays(end).name=[params.EVnames{i},'_fdr_Z (radius = ',num2str(params.radius),')'];
            for iScan = params.scanNum
                overlays(end).data{iScan}=out{iScan}.fdrDiag_Class(:,:,:,i);
                overlays(end).params=params.scanParams{iScan};
            end
        end
    end
    if params.fweAdjustment
        overlays(end+1)=thisOverlay;
        overlays(end).name=['fwe Z (radius = ',num2str(params.radius),')'];
        for iScan = params.scanNum
            overlays(end).data{iScan}=out{iScan}.fweDiag;
            overlays(end).params=params.scanParams{iScan};
        end
        for i=1:params.numberEvents
            overlays(end+1)=thisOverlay;
            overlays(end).name=[params.EVnames{i},'_fwe_Z (radius = ',num2str(params.radius),')'];
            for iScan = params.scanNum
                overlays(end).data{iScan}=out{iScan}.fweDiag_Class(:,:,:,i);
                overlays(end).params=params.scanParams{iScan};
            end
        end
    end
end

overlays(end+1)=defaultOverlay;
overlays(end).name = ['searchlight size (radius = ',num2str(params.radius),')'];
maxValue = 0;
minValue = inf;
for iScan = params.scanNum
  maxValue = max(maxValue,nanmax(out{iScan}.size_light(:)));
  minValue = min(minValue,nanmin(out{iScan}.size_light(:)));
  overlays(end).data{iScan}=out{iScan}.size_light;
  overlays(end).params{iScan} = params.scanParams{iScan};
end
overlays(end).clip = [minValue maxValue];
overlays(end).range =overlays(end).clip;

classAnal.overlays = overlays;
thisView = viewSet(thisView,'newAnalysis',classAnal);
if ~isempty(viewGet(thisView,'fignum'))
  refreshMLRDisplay(viewGet(thisView,'viewNum'));
end

%-------------------------------------------------------- Save the analysis
saveAnalysis(thisView,classAnal.name);

set(viewGet(thisView,'figNum'),'Pointer','arrow');

return
