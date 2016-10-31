% scanInfo.m
%
%      usage: scanInfo(scanNum,groupNum,<displayInDialog>)
%         by: justin gardner
%       date: 11/08/06
%    purpose: print out information about scans in group
%       e.g.: scanInfo(1);
%             scanInfo(3,'Raw');
function retval = scanInfo(scanNum,groupNum,displayInDialog)

% check arguments
if ~any(nargin == [2 3])
  help scanInfo
  return
end

if ~exist('displayInDialog','var'),displayInDialog = 0;end

view = newView;

% if groupNum is a string, then user passed in a name rather than
% a number
if isstr(groupNum)
  groupName = groupNum;
  groupNum = viewGet(view,'groupNum',groupName);
  if isempty(groupNum)
    disp(sprintf('(groupInfo): No group %s',groupName));
    deleteView(view);
    return
  end
end

if (groupNum < 1) | (groupNum > viewGet(view,'numberOfGroups'))
  disp(sprintf('(groupInfo): No group number %i',groupNum));
  deleteView(view);
  return
end
groupName = viewGet(view,'groupName',groupNum);
% set the group and scan
view = viewSet(view,'curGroup',groupNum);

if (scanNum < 1) || (scanNum > viewGet(view,'nScans'))
  disp(sprintf('(scanInfo) Could not find scan %i in group %s',scanNum,groupName));
  deleteView(view);
  return
end

% grab info
disppercent(-inf,'(scanInfo) Gathering scan info');
description = viewGet(view,'description',scanNum,groupNum);
scanVoxelSize = viewGet(view,'scanVoxelSize',scanNum,groupNum);
tr = viewGet(view,'TR',scanNum,groupNum);
framePeriod = viewGet(view,'framePeriod',scanNum,groupNum);
totalFrames = viewGet(view,'totalFrames',scanNum,groupNum);
filename = viewGet(view,'tSeriesFile',scanNum,groupNum);
originalFilename = viewGet(view,'originalFilename',scanNum,groupNum);
originalGroupname = viewGet(view,'originalGroupname',scanNum,groupNum);
stimFilename = cellArray(viewGet(view,'stimFilename',scanNum,groupNum),-1);
scanHdr = viewGet(view,'niftiHdr',scanNum,groupNum);
scanQform = viewGet(view,'scanQform',scanNum,groupNum);
scanQformCode = viewGet(view,'scanQformCode',scanNum,groupNum);
scanSform = viewGet(view,'scanSform',scanNum,groupNum);
scanSformCode = viewGet(view,'scanSformCode',scanNum,groupNum);
scanDims = viewGet(view,'scanDims',scanNum,groupNum);
junkFrames = viewGet(view,'junkFrames',scanNum,groupNum);
totalJunkedFrames = viewGet(view,'totalJunkedFrames',scanNum,groupNum);
vol2mag = viewGet(view,'scanVol2mag',scanNum,groupNum);
vol2tal = viewGet(view,'scanVol2tal',scanNum,groupNum);

matfile = viewGet(view,'params',scanNum,groupNum);
fieldsToIgnore = {'tseriesfiles','groupname','description','filename'};
extraFields = {};extraFieldsValue = {};
% get extra parameters form associated matfile
if ~isempty(matfile) && isfield(matfile,'params')
  fields = fieldnames(matfile.params);
  for i = 1:length(fields)
    % ignore fields from the list set above
    if ~any(strcmp(lower(fields{i}),fieldsToIgnore))
      % and fields that have groupname in them
      if isempty(strfind(lower(fields{i}),'groupname'))
	% display fields that are numeric
	if isnumeric(matfile.params.(fields{i}))
	  extraFields{end+1} = fields{i};
	  extraFieldsValue{end+1} = num2str(matfile.params.(fields{i}));
	% or are strings
	elseif isstr(matfile.params.(fields{i}))
	  extraFields{end+1} = fields{i};
	  extraFieldsValue{end+1} = matfile.params.(fields{i});
	end
      end
    end
  end
end

% get params info
nCols = 1;
% display info
if displayInDialog | (nargout == 1)
  paramsInfo = {{'description',description,'editable=0','Scan description'},...
		{'Filename',filename,'editable=0','Name of file'},...
		{'GroupName',groupName,'editable=0','Name of group'}};
  for i = 1:length(originalFilename)
    paramsInfo{end+1} = {sprintf('Original%i',i) sprintf('%s: %s',originalGroupname{i},originalFilename{i}) 'editable=0' 'Name of original group and filename that this scan came from'};
  end
  if ~isempty(stimFilename)
    % make stimfiles into sets of 4 (otherwise it takes too much room in the
    % dialog for long concatenations or averages).
    for j = 1:ceil(length(stimFilename)/4)
      stimfileNames = getLastDir(stimFilename{(j-1)*4+1});
      firstStimfile = ((j-1)*4+1);
      lastStimfile = min(length(stimFilename),j*4);
      for i = firstStimfile+1:lastStimfile
	stimfileNames = sprintf('%s, %s',stimfileNames,getLastDir(stimFilename{i}));
      end
      if lastStimfile == firstStimfile
	paramsInfo{end+1} = {sprintf('stimFilenames%i',lastStimfile),stimfileNames,'editable=0','Names of associated stimfiles'};
      else
	paramsInfo{end+1} = {sprintf('stimFilenames%ito%i',firstStimfile,lastStimfile),stimfileNames,'editable=0','Names of associated stimfiles'};
      end
    end
  end
  paramsInfo{end+1} = {'voxelSize',scanVoxelSize,'editable=0','Voxel dimensions in mm'};
  paramsInfo{end+1} =  {'dims',scanDims,'editable=0','Dimensions of scan'};
  paramsInfo{end+1} =  {'TR',tr,'editable=0','TR. This is retrieved frome the dicom header information.'};
  paramsInfo{end+1} =  {'framePeriod',framePeriod,'editable=0','Time each volume takes to acquire. Note that this is usually the same as TR, except for 3D scans in which this should be number of slices times the TR.'};
  paramsInfo{end+1} =  {'numVolumes',totalFrames,'editable=0','Number of volumes'};
  paramsInfo{end+1} =  {'junkFrames',junkFrames,'editable=0','Number of junk frames'};
  paramsInfo{end+1} =  {'totalJunkedFrames',num2str(totalJunkedFrames),'type=string','editable=0','Number of junk frames that have already been discarded from this time series (this is useful for aligning the volume number of a stimulus set in a stimulus file with the volumes in the time series)'};
  paramsInfo{end+1} = {'qform',scanQform,'editable=0','Qform matrix specifies the transformation to the scanner coordinate frame'};
  paramsInfo{end+1} = {'qformCode',scanQformCode,'editable=0','If qformCode is 0 it means the qform has never been set.'};
  paramsInfo{end+1} = {'sform',scanSform,'editable=0','Sform matrix is set by mrAlign and usually specifies the transformation to the volume anatomy'};
  paramsInfo{end+1} = {'sformCode',scanSformCode,'editable=0','If sformCode is 0 it means the sform has never been set and mrLoadRet will use the qform to compute the transform to the base anatomy. If mrAlign has been run properly, then this value should be set to 1. This value can also be 3 if the sform is a transformation to talairach coordinates.'};
  paramsInfo{end+1} = {'vol2mag',vol2mag,'editable=0','This is the xformation that takes the canonical base that this scan was aligned to into magnet coordinates. To set this field, you can use export from mrAlign.'};
  paramsInfo{end+1} = {'vol2tal',vol2tal,'editable=0','This is the xformation that takes the canonical base that this scan was aligned to into talairach coordinates. To set this field, you can export talairach coordinates from mrAlign.'};

  % display extra parameters form associated matfile
  for i = 1:length(extraFields)
    paramsInfo{end+1} = {extraFields{i} extraFieldsValue{i} 'editable=0','Parameter from associated mat file'};
  end

  disppercent(inf);
  if (nargout == 1)
    retval = mrParamsDefault(paramsInfo);
  else
    mrParamsDialog(paramsInfo,'Scan info',nCols);
  end
else
  disppercent(inf);
  disp(sprintf('%s',description));
  disp(sprintf('Filename: %s GroupName: %s',filename,groupName));
  for i = 1:length(originalFilename)
    disp(sprintf('Original Filename: %s Group: %s',originalFilename{i},originalGroupname{i}));
  end
  for i = 1:length(stimFilename)
    disp(sprintf('StimFilename: %s',stimFilename{i}));
  end

  if isempty(tr),tr=nan;end
  disp(sprintf('voxelSize=[%0.1f %0.1f %0.1f] TR=%0.4f framePeriod=%0.4f Dims: [%i %i %i] Volumes=%i',scanVoxelSize(1),scanVoxelSize(2),scanVoxelSize(3),tr,framePeriod,scanDims(1),scanDims(2),scanDims(3),totalFrames));
  disp(sprintf('junkFrames=%i totalJunkedFrames=[%s]',junkFrames,num2str(totalJunkedFrames)));

  % display qform and sform
  disp(sprintf('++++++++++++++++++++++++++ qform ++++++++++++++++++++++++++'));
  for rownum = 1:4
    disp(sprintf('%f\t%f\t%f\t%f',scanQform(rownum,1),scanQform(rownum,2),scanQform(rownum,3),scanQform(rownum,4)));
  end
  disp(sprintf('qformCode: %i',scanQformCode));

  disp(sprintf('++++++++++++++++++++++++++ sform ++++++++++++++++++++++++++'));
  for rownum = 1:4
    disp(sprintf('%f\t%f\t%f\t%f',scanSform(rownum,1),scanSform(rownum,2),scanSform(rownum,3),scanSform(rownum,4)));
  end
  
  % display extra fields
  for i = 1:length(extraFields)
    disp(sprintf('%s: %s',extraFields{i},extraFieldsValue{i}));
  end
  
  disp(sprintf('sformCode: %i',scanSformCode));
end

deleteView(view);


