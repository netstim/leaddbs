% eventRelated.m
%
%      usage: view = eventRelated(view,params)
%         by: justin gardner
%       date: 10/20/06
%    purpose: event related data analysis
%
%             if you just want a default parameter structure you
%             can do:
% 
%             v = newView;
%             [v params] = eventRelated(v,[],'justGetParams=1','defaultParams=1','scanList=1')
%
%             Note that justGetParams,defualtParams and scanList are independent parameters, so
%             if you want, say to bring up the GUI to set the params, but not run the analysis, you
%             can do:
%             [v params] = eventRelated(v,[],'justGetParams=1');
%
function [view d] = eventRelated(view,params,varargin)

d = [];

% check arguments
if ~any(nargin == [1 2 3 4 5 6 7 8])
  help eventRelated.m
  return
end

mrGlobals;

% other arguments
eval(evalargs(varargin,[],[],{'justGetParams','defaultParams','scanList'}));
if ieNotDefined('justGetParams'),justGetParams = 0;end
if ieNotDefined('defaultParams'),defaultParams = 0;end
if ieNotDefined('scanList'),scanList = [];end

% First get parameters
if ieNotDefined('params')
  % put up the gui
  if defaultParams
    params = eventRelatedGUI('groupName',viewGet(view,'groupName'),'useDefault=1','scanList',scanList);
  else
    params = eventRelatedGUI('groupName',viewGet(view,'groupName'),'scanList',scanList);
  end
end

% check params (and possibly reformat to latest version)
params = checkEventRelatedParams(params);
  
% just return parameters
if justGetParams
  d = params;
  return
end

% Reconcile params with current status of group and ensure that it has
% the required fields. 
params = mrParamsReconcile([],params);

% Abort if params empty
if ieNotDefined('params'),return,end

% set the group
view = viewSet(view,'groupName',params.groupName);

if ~isfield(params, 'applyFiltering')
  params.applyFiltering = 0;
end

% inplace concatenation is handeled by a different function
if isfield(params, 'inplaceConcat')
    if params.inplaceConcat
        [view d] = eventRelatedMultiple(view,params);
        return
    end
end

% create the parameters for the overlay
dateString = datestr(now);
r2.name = 'r2';
r2.groupName = params.groupName;
r2.function = 'eventRelated';
r2.reconcileFunction = 'defaultReconcileParams';
r2.data = cell(1,viewGet(view,'nScans'));
r2.date = dateString;
r2.params = cell(1,viewGet(view,'nScans'));
r2.range = [0 1];
r2.clip = [0 1];
% colormap is made with a little bit less on the dark end
r2.colormap = hot(312);
r2.colormap = r2.colormap(end-255:end,:);
r2.alpha = 1;
r2.colormapType = 'setRangeToMax';
r2.interrogator = 'eventRelatedPlot';
r2.mergeFunction = 'defaultMergeParams';

tic
set(viewGet(view,'figNum'),'Pointer','watch');drawnow;
for scanNum = params.scanNum
  % decide how many slices to do at a time, this is done
  % simply to save memory -- currently our system is limited
  % to 2G of memory and for large concatenations, you need
  % to break up the analysis into smaller portions of the data
  numSlices = viewGet(view,'nSlices',scanNum);
  numVolumes = viewGet(view,'nFrames',scanNum);
  dims = viewGet(view,'dims',scanNum);
  % choose how many slices based on trying to keep a certain
  % amount of data in the memory
  [numSlicesAtATime rawNumSlices numRowsAtATime precision] = getNumSlicesAtATime(numVolumes,dims);
  currentSlice = 1;
  ehdr = [];ehdrste = [];thisr2 = [];

  for i = 1:ceil(numSlices/numSlicesAtATime)
    % calculate which slices we are working on
    thisSlices = [currentSlice min(numSlices,currentSlice+numSlicesAtATime-1)];
    % set the row we are working on
    currentRow = 1;
    % clear variables that will hold the output for the slices we
    % are working on
    sliceEhdr = [];sliceEhdrste = [];sliceR2 = [];
    for j = 1:ceil(dims(1)/numRowsAtATime)
      % load the scan
      thisRows = [currentRow min(dims(1),currentRow+numRowsAtATime-1)];
      d = loadScan(view,scanNum,[],thisSlices,precision,thisRows);
      % get the stim volumes, if empty then abort
      d = getStimvol(d,params.scanParams{scanNum}.stimvolVarInfo);
      if isempty(d.stimvol),mrWarnDlg('No stim volumes found');return,end
      % do any called for preprocessing
      d = eventRelatedPreProcess(d,params.scanParams{scanNum}.preprocess);
      % make a stimulation convolution matrix
      d = makescm(d,ceil(params.scanParams{scanNum}.hdrlen/d.tr),params.applyFiltering);
      % compute the estimated hemodynamic responses
      d = getr2(d);
      % update the current row we are working on
      currentRow = currentRow+numRowsAtATime;
      % if we are calculating full slice, then just pass that on
      if numRowsAtATime == dims(1)
	sliceEhdr = d.ehdr;
	sliceEhdrste = d.ehdrste;
	sliceR2 = d.r2;
      % working on a subset of rows, cat together with what
      % has been computed for other rows
      else
	sliceEhdr = cat(1,sliceEhdr,d.ehdr);
	sliceEhdrste = cat(1,sliceEhdrste,d.ehdrste);
	sliceR2 = cat(1,sliceR2,d.r2);
      end
    end
    % update the current slice we are working on
    currentSlice = currentSlice+numSlicesAtATime;
    % cat with what has already been computed for other slices
    ehdr = cat(3,ehdr,sliceEhdr);
    ehdrste = cat(3,ehdrste,sliceEhdrste);
    thisr2 = cat(3,thisr2,sliceR2);
  end

  % now put all the data from all the slices into the structure
  d.ehdr = single(ehdr);
  d.ehdrste = single(ehdrste);
  d.r2 = single(thisr2);

  % get the actual size of the data (not just the size of the last
  % slice/set of rows we were working on).
  d.dim(1:3) = size(d.r2);

  % save the r2 overlay
  r2.data{scanNum} = d.r2;
  r2.params{scanNum} = params.scanParams{scanNum};
  
  % save other eventRelated parameters
  erAnal.d{scanNum}.ver = d.ver;
  erAnal.d{scanNum}.filename = d.filename;
  erAnal.d{scanNum}.filepath = d.filepath;
  erAnal.d{scanNum}.dim = d.dim;
  erAnal.d{scanNum}.ehdr = d.ehdr;
  erAnal.d{scanNum}.ehdrste = d.ehdrste;
  erAnal.d{scanNum}.nhdr = d.nhdr;
  erAnal.d{scanNum}.hdrlen = d.hdrlen;
  erAnal.d{scanNum}.tr = d.tr;
  erAnal.d{scanNum}.stimvol = d.stimvol;
  erAnal.d{scanNum}.stimNames = d.stimNames;
  erAnal.d{scanNum}.scm = d.scm;
  erAnal.d{scanNum}.expname = d.expname;
  erAnal.d{scanNum}.fullpath = d.fullpath;
end
toc

% install analysis
erAnal.name = params.saveName;
erAnal.type = 'erAnal';
erAnal.groupName = params.groupName;
erAnal.function = 'eventRelated';
erAnal.reconcileFunction = 'defaultReconcileParams';
erAnal.mergeFunction = 'defaultMergeParams';
erAnal.guiFunction = 'eventRelatedGUI';
erAnal.params = params;
erAnal.overlays = r2;
erAnal.curOverlay = 1;
erAnal.date = dateString;
view = viewSet(view,'newAnalysis',erAnal);
if ~isempty(viewGet(view,'fignum'))
  refreshMLRDisplay(viewGet(view,'viewNum'));
end

% Save it
saveAnalysis(view,erAnal.name);

set(viewGet(view,'figNum'),'Pointer','arrow');drawnow

% for output
if nargout > 1
  for i = 1:length(d)
    erAnal.d{i}.r2 = r2.data{i};
  end
  % make d strucutre
  if length(erAnal.d) == 1
    d = erAnal.d{1};
  else
    d = erAnal.d;
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    checkEventRelatedParams    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function params = checkEventRelatedParams(params)
  
if isempty(params)
  disp(sprintf('(eventRelated:checkEventRelatedParams) Empty params'));
  return
end

% convert scanParams if necessary from old to new style. New style
% has one field call stimvolVarInfo which contains the varname, taskNum, segmentNum and phaseNum
% for computing stimvols using getStimvol. This was done so that stimvolVarInfo can
% be a cell array with multiple variable settings in it (which will get concatenated)
for scanNum = params.scanNum
  if ~isempty(params.scanParams{scanNum}) && ~isfield(params.scanParams{scanNum},'stimvolVarInfo')
    fieldsToMove = {'phaseNum','taskNum','segmentNum','varname','stimtrace'};
    for iField = 1:length(fieldsToMove)
      if isfield(params.scanParams{scanNum},fieldsToMove{iField})
	% copy into the combined stimvolVarInfo field
	params.scanParams{scanNum}.stimvolVarInfo.(fieldsToMove{iField}) = params.scanParams{scanNum}.(fieldsToMove{iField});
	% and remove from original location
	params.scanParams{scanNum} = rmfield(params.scanParams{scanNum},fieldsToMove{iField});
      end
    end
  end
  % check for old style
  if ~isfield(params.scanParams{scanNum},'stimvolVarInfo')
    % this should not happen - it means that none of the usual fields for
    % what stimulus event to trigger off of have been set.
    disp(sprintf('(eventRelated:checkEventRelatedParams) No mgl stimfile information available.'));
    params.scanParams{scanNum}.stimvolVarInfo = [];
  end
end