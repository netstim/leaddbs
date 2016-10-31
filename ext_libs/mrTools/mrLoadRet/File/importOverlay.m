% importOverlay.m
%
%        $Id$ 
%      usage: importOverlay(thisView)
%         by: julien besle
%       date: 20/01/2010
%    purpose: import overlay in an analysis of mrLoadRet GUI
%
function retval = importOverlay(thisView)

% check arguments
if ~any(nargin == [1])
  help importOverlay
  return
end

if ~isanalysis(viewGet(thisView,'analysis'))
    mrWarnDlg('(importOverlay) Overlays must be imported into an analysis. Use Edit -> Analysis -> New Analysis.')
    return
end

%first check that there is at least one scan and display scan selection dialog if more than one
scan_number = viewGet(thisView,'nScans');
switch(scan_number)
  case 0
    mrWarnDlg('(importOverlay) Group must have at least one scan.')
    return
  case 1
    scanlist = 1;
  otherwise
    scanlist = selectInList(thisView,'scans');
end

% go find the group that user wants to load here
[filename pathname] = uigetfile({sprintf('*%s',mrGetPref('niftiFileExtension')),'Nifti files'},'Select nifti file that you want to import');
if (filename==0)
  return
end

% get the full file name
fullFilename = fullfile(pathname,filename);

% read the nifti header
[data,hdr] = mlrImageReadNifti(fullFilename);
if isempty(hdr),return,end

% make sure it has only 1 frame
if hdr.dim(1)==3
   hdr.dim(5)=1;
end
if hdr.dim(5) < 1
  mrWarnDlg(sprintf('(importOverlay) Could not import image because it has %d frames',hdr.dim(5)));
  return
elseif hdr.dim(5) > 1
  frameList = buttondlg('Choose the frames to import', [cellstr(int2str((1:hdr.dim(5))'));{'All'}]);
  if isempty(frameList)
    return
  end
  if frameList(end)
    frameList = 1:hdr.dim(5);
  else
    frameList = find(frameList);
  end
else
  frameList=1;
end
nFrames = length(frameList);



for i_scan = scanlist
scan_dims = viewGet(thisView,'datasize',i_scan);
   if any(hdr.dim([2 3 4])'~= scan_dims)
      mrWarnDlg(sprintf('(importOverlay) Could not import image because it is not compatible with scan (%s vs %s)', num2str(hdr.dim([2 3 4])'),num2str(scan_dims)) );
      return
   end
end
   
maxFrames=20;
paramsInfo = {{'filename',filename,'The name of the nifti file that you are importing','editable=0'}};
if nFrames<=maxFrames
  for iFrame = 1:nFrames
    numFrame = num2str(frameList(iFrame));
    paramsInfo{end+1} = {['nameFrame' numFrame],['Frame ' numFrame],['A description for the overlay corresponding to frame' numFrame]};
  end
else
  paramsInfo{end+1} = {'nameFrame','','A description for the overlay you are importing. the frame number will be added at the end for each overlay'};
end
params = mrParamsDialog(paramsInfo);
if isempty(params),return,end
if nFrames>maxFrames
  for iFrame = 1:nFrames
    numFrame = num2str(frameList(iFrame));
    params.(['nameFrame' numFrame]) = [params.nameFrame '(Frame ' numFrame ')'];
  end
end
  
max_overlay = max(data(:));
min_overlay = min(data(:));
defaultOverlay.name = '';
defaultOverlay.groupName = viewGet(thisView,'groupName');
defaultOverlay.function = 'importOverlay';
defaultOverlay.reconcileFunction = 'defaultReconcileParams';
defaultOverlay.data = cell(scan_number,1);
defaultOverlay.date = datestr(now);
defaultOverlay.params = params;
% colormap is made with a little bit less on the dark end
defaultOverlay.colormap = hot(312);
defaultOverlay.colormap = defaultOverlay.colormap(end-255:end,:);
defaultOverlay.alpha = 1;
defaultOverlay.interrogator = '';
defaultOverlay.mergeFunction = 'defaultMergeParams';
defaultOverlay.colormapType = 'normal';
defaultOverlay.range = [min_overlay max_overlay];
defaultOverlay.clip = [min_overlay max_overlay];

for iFrame=1:nFrames
  numFrame = frameList(iFrame);
  overlays(iFrame) = defaultOverlay;
  overlays(iFrame).name = [params.(['nameFrame' num2str(numFrame)]) ' (' params.filename ')']; 
  overlays(iFrame).data(scanlist) = squeeze(num2cell(repmat(data(:,:,:,numFrame),[1 1 1 length(scanlist)]),[1 2 3]))';
end

thisView = viewSet(thisView,'newoverlay',overlays);
refreshMLRDisplay(thisView.viewNum);





