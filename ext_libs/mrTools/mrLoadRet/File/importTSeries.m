% importTSeries.m
%
%        $Id$ 
%      usage: importTSeries(v)
%         by: justin gardner
%       date: 07/03/09
%    purpose: 
%
function retval = importTSeries(v)

% check arguments
if ~any(nargin == [1])
  help importTSeries
  return
end

minFramePeriod = .01;  %frame period in sec outside which the user is prompted
maxFramePeriod = 100;  % that something weird's goin on

% go find the group that user wants to load here
[filename pathname] = uigetfile({'*.nii;*.img','Nifti files'},'Select nifti tSeries that you want to import','multiselect','on');

if isnumeric(filename)
  return
elseif ~iscell(filename)
  filename = {filename};
end

nScans = viewGet(v,'nScans');
cScan=0;  %count of succesfully added scans
weirdFramePeriods = 0;

for iFile = 1:length(filename)
  
  if ~isempty(strfind(stripext(filename{iFile}),'.'))
    [~,name,extension] = fileparts(filename{iFile});
    mrWarnDlg(sprintf('(importTSeries) Ignoring file %s because it has a . in the filename that does not mark the file extension. If you want to use this file, consider renaming to %s',filename{iFile},setext(fixBadChars(name,{'.','_'}),extension)));
    return
  end

  % get the full file name
  fullFilename = fullfile(pathname,filename{iFile});

  % read the nifti header
  hdr = mlrImageReadNiftiHeader(fullFilename);
  if isempty(hdr),return,end
  % make sure it is 4 dimensional and then get the number of frames
  if ~ismember(hdr.dim(1),[3 4])
    mrWarnDlg(sprintf('(importTSeries) Could not import tSeries because it is a %i dimensional file',hdr.dim(1)));
    return
  elseif hdr.dim(1)==3
    nFrames = 1;
  else
    nFrames = hdr.dim(5);
  end

  paramsInfo = {{'filename',filename{iFile},'The name of the nifti file that you are importing','editable=0'}};
  paramsInfo{end+1} = {'description','','A description for the nifti tSeries you are imporint'};
  paramsInfo{end+1} = {'nFrames',nFrames,'incdec=[-1 1]',sprintf('minmax=[0 %i]',nFrames),'Number of total frames in your nfiti tSeries'};
  paramsInfo{end+1} = {'junkFrames',0,'incdec=[-1 1]',sprintf('minmax=[0 %i]',nFrames),'How many frames should be junked at the beginning'};

  params = mrParamsDialog(paramsInfo);
  if isempty(params),return,end

  thisScanParams.fileName = params.filename;
  thisScanParams.description = params.description;
  thisScanParams.nFrames = params.nFrames;
  thisScanParams.junkFrames = params.junkFrames;

  % now read the file
  %tSeries = mlrImageReadNifti(fullFilename);

  v = saveNewTSeries(v,fullFilename,thisScanParams);
  
  %get the new scan params and check that frame period is a reasonable value
  newScanNum = viewGet(v,'nScans');
  if newScanNum>nScans %if the new scan has been succesfully added
    cScan=cScan+1;
    scanNum(cScan)=newScanNum;
    scanParams(cScan) = viewGet(v,'scanparams',newScanNum);
    if scanParams(cScan).nFrames>1 && (scanParams(cScan).framePeriod>maxFramePeriod || scanParams(cScan).framePeriod<minFramePeriod)
      weirdFramePeriods(cScan) = scanParams(cScan).framePeriod;
    else
      weirdFramePeriods(cScan) = 0;
    end
  end
  
end

%fix apparently abnormal frame periods
if any(weirdFramePeriods)
  scanParams = fixFramePeriods(scanParams,weirdFramePeriods,minFramePeriod,maxFramePeriod,viewGet(v,'tseriesDir'));
end

%set the new frame period values
for iScan = 1:length(scanNum)
  if weirdFramePeriods(iScan)
    v=viewSet(v,'scanParams',scanParams(iScan),scanNum(iScan));
  end
end



