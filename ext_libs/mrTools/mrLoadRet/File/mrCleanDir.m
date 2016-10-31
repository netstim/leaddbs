% mrCleanDir.m
%
%        $Id$
%      usage: mrCleanDir(<force>)
%         by: justin gardner
%       date: 10/20/06
%    purpose: Removes (or recovers) scans that have been deleted
%             using deleteScans. DeleteScans simpy removes the scan
%             from the MLR mrSession variable, but does not destroy
%             the original file. This deletes the nifti files
%             Call with force=1 if you want to remove all files
%             without being asked what to do. Default is to be asked
% 
%
function retval = mrCleanDir(force)

% check arguments
if ~any(nargin == [0 1])
  help mrCleanDir
  return
end

if ieNotDefined('force'),force = 0;end

view = newView;

% check for unlinked files
groups = viewGet(view,'groupNames');

for g = 1:length(groups)
  % set the current group
  groupNum = viewGet(view,'groupNum',groups{g});
  view = viewSet(view,'curGroup',groupNum);
  nScans = viewGet(view,'nScans');
  
  if nScans > 0
    % get the directory
    [tseriesDirName] = fileparts(viewGet(view,'tseriesPath',1));
  else
    tseriesDirName = fullfile(viewGet(view,'groupName'),'TSeries');
  end
  tseriesDir = [dir(sprintf('%s/*.hdr',tseriesDirName));dir(sprintf('%s/*.nii',tseriesDirName))];
  for i = 1:length(tseriesDir)
    tseriesDir(i).match = 0;
  end

  
  % look for unmatched files
  for scanNum = 1:nScans
    tseriesFilename = viewGet(view,'tseriesPath',scanNum,groupNum);
    [thisDirName thisFilename] = fileparts(tseriesFilename);
    for filenum = 1:length(tseriesDir)
      [dirname,filename] = fileparts(tseriesDir(filenum).name);
      if strcmp(filename,thisFilename) && strcmp(tseriesDirName,thisDirName)
	tseriesDir(filenum).match = 1;
      end
    end
  end

  % count to see if we have all matches
  matched = 0;
  for i = 1:length(tseriesDir)
    if tseriesDir(i).match
      matched = matched+1;
    end
  end

  % if we have more files in directory than that are matched,...
  if (matched < length(tseriesDir))
    disp(['%%%%%%%%%%%%%%% Group ',groups{g}]);
    recoverable = [];
    % display the names of the hdr/img/mat files that are not matched
    for i = 1:length(tseriesDir)
      disp(sprintf('================ScanNum %i =============================',i));
      [path baseFilename] = fileparts(tseriesDir(i).name);
      [recoverable(i) scanParams{i}] = dispParams(tseriesDirName,baseFilename,tseriesDir(i).match);
      if tseriesDir(i).match
	disp(sprintf('Matched'));
      end
    end
    % see if there are any recoverable files
    if sum(recoverable)
      for i = 1:length(tseriesDir)
	if recoverable(i) && ~force
	  if askuser(sprintf('Recover scan %i',i))
	    disp(sprintf('Recovering scan %i',i));
	    view = viewSet(view,'newScan',scanParams{i});
	    tseriesDir(i).match = 1;
	  end
	end
      end
    end
    % and ask user if they should be deleted
    if force | askuser(sprintf('Delete files from group %s',groups{g}))
      for i = 1:length(tseriesDir)
	if ~tseriesDir(i).match
	  [path baseFilename] = fileparts(tseriesDir(i).name);
	  filename = sprintf('%s/%s.hdr',tseriesDirName,baseFilename);
	  if isfile(filename),delete(filename),end;
	  filename = sprintf('%s/%s.img',tseriesDirName,baseFilename);
	  if isfile(filename),delete(filename),end;
	  filename = sprintf('%s/%s.nii',tseriesDirName,baseFilename);
	  if isfile(filename),delete(filename),end;
	  filename = sprintf('%s/%s.mat',tseriesDirName,baseFilename);
	  if isfile(filename),delete(filename),end;
	  
  end
      end
    end
  elseif (nScans > length(tseriesDir))
    disp(['!!!! Missing scans in group ' groups{g} ' (' num2str(length(tseriesDir)) ':' num2str(nScans) ') !!!!']);
  else
    disp(sprintf('Group %s matches (%i:%i)',groups{g},length(tseriesDir),nScans));
  end

end

deleteView(view);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% display the parameters and filename for the unlinked scan
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [recoverable scanParams] = dispParams(tseriesDirName,baseFilename,match)

recoverable = 0;scanParams = [];
% see if there is a mat file
filename = fullfile(tseriesDirName,sprintf('%s.mat',baseFilename));
if isfile(filename)
  matfile = load(filename);
  % first see if we can match the tseriesFileName with
  % the one we have here
  if isfield(matfile,'params')
    if isfield(matfile,'tseriesFileName') 
      tseriesFileName = fullfile(tseriesDirName,matfile.tseriesFileName);
      tseriesHdrFileName = sprintf('%s.hdr',stripext(tseriesFileName));
    else
      tseriesFileName = sprintf('%s.img',fullfile(tseriesDirName,baseFilename));
      tseriesHdrFileName = sprintf('%s.hdr',fullfile(tseriesDirName,baseFilename));
    end
    % check if they are there
    if isfile(tseriesFileName) && isfile(tseriesHdrFileName)
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % if this is a motionComp then display its parameters
      % and see if we can recover it
      if isfield(matfile.params,'motionCompGroupName')
	persistent motionCompParams;
	% display the parameters (only if this is one that hasn't
	% been matched, otherwise just keep a record of that
	% param file.
	if ~match
	  dispMotionCompParams(matfile.params);
	  % check how many times we have seen this
	  targetNum = sum(structIsMember(matfile.params,motionCompParams))+1;
	  if length(matfile.params.descriptions) >= targetNum
	    disp(sprintf('Description: %s',matfile.params.descriptions{targetNum}));
	    % get a view and set it to the original group
	    v = newView;
	    v = viewSet(v,'curGroup',viewGet(v,'groupNum',matfile.params.groupName));
	    % read the image header
	    hdr = mlrImageReadNiftiHeader(tseriesHdrFileName);
	    % get the scan params
	    scanParams.junkFrames = 0;
	    scanParams.nFrames = hdr.dim(5);
	    scanParams.description = matfile.params.descriptions{targetNum};
	    scanParams.fileName = getLastDir(tseriesFileName);
	    scanParams.originalFileName{1} = viewGet(v,'tSeriesFile',matfile.params.targetScans(targetNum));
	    scanParams.originalGroupName{1} = matfile.params.groupName;
	    %set that we can recover this file
	    recoverable = 1;
	  end
	end
	motionCompParams{end+1} = matfile.params;
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      else
	if ~match
	  % display as average
	  if isfield(matfile.params,'aveGroupName')
	    dispAverageParams(matfile.params);
	  % or as a set of files
	  else 
	    filename = sprintf('%s/%s.hdr',tseriesDirName,baseFilename);
	    if isfile(filename),disp(filename),end
	    filename = sprintf('%s/%s.img',tseriesDirName,baseFilename);
	    if isfile(filename),disp(filename),end
	    filename = sprintf('%s/%s.nii',tseriesDirName,baseFilename);
	    if isfile(filename),disp(filename),end
	    filename = sprintf('%s/%s.mat',tseriesDirName,baseFilename);
	    if isfile(filename),disp(filename),end
	  end
	  % read the image header
	  hdr = mlrImageReadNiftiHeader(tseriesHdrFileName);
	  scanParams.junkFrames = 0;
	  scanParams.nFrames = hdr.dim(5);
	  if isfield(matfile.params,'description')
	    disp(sprintf('description: %s',matfile.params.description));
	    scanParams.description = matfile.params.description;
	  else
	    scanParams.description = '';
	  end	    
	  scanParams.fileName = getLastDir(tseriesFileName);
	  % get original file and group
	  if isfield(matfile.params,'tseriesfiles')
	    originalFiles = matfile.params.tseriesfiles;
	  elseif isfield(matfile.params,'tseriesFile')
	    originalFiles = matfile.params.tseriesFile;
	  else
	    originalFiles = {};
	  end
	  if isfield(matfile.params,'groupName')
	    originalGroup = matfile.params.groupName;
	  else
	    originalGroup = 'Unknown';
	  end
	  for scanNum = 1:length(originalFiles)
	    scanParams.originalFileName{scanNum} = originalFiles{scanNum};
	    scanParams.originalGroupName{scanNum} = originalGroup;
	  end
	  recoverable = 1;
	end
      end

    end
  end
end

% if we don't know how to recover then just show files for deleting
if ~recoverable && ~match
  filename = sprintf('%s/%s.hdr',tseriesDirName,baseFilename);
  if isfile(filename),disp(filename),end
  filename = sprintf('%s/%s.img',tseriesDirName,baseFilename);
  if isfile(filename),disp(filename),end
  filename = sprintf('%s/%s.nii',tseriesDirName,baseFilename);
  if isfile(filename),disp(filename),end
  filename = sprintf('%s/%s.mat',tseriesDirName,baseFilename);
  if isfile(filename),disp(filename),end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% display the parameters in a motion comp params structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dispMotionCompParams(params)

if isfield('params', 'correctIntensityContrast') 
  % this must be an old motionComp params.  that field doesn't exist anymore
  disp(sprintf('GroupName: %s baseScan: %i baseFrame: %s robust: %i correctIntensityContrast: %i',params.groupName,params.baseScan,params.baseFrame,params.robust,params.correctIntensityContrast));
else
  disp(sprintf('GroupName: %s baseScan: %i baseFrame: %s robust: %i gradIntensityCorrection: %i',params.groupName,params.baseScan,params.baseFrame,params.robust,params.gradIntensityCorrection));
end

disp(sprintf('crop: %s niters: %i interpMethod: %s targetScans: %s',num2str(params.crop),params.niters,params.interpMethod,num2str(params.targetScans)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% display the parameters in a motion comp params structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dispAverageParams(params)

disp(sprintf('Average from GroupName: %s baseScan: %i',params.groupName,params.baseScan));
for i = 1:length(params.tseriesfiles)
  disp(sprintf('%i %s: %s Shift: %i Reverse: %i',params.scanList(i),params.groupName,params.tseriesfiles{i},params.shiftList(i),params.reverseList(i)));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find the matching struct
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function retval = structIsMember(s,list)

retval = [];
for i = 1:length(list)
  if isequal(s,list{i})
    retval(i) = 1;
  else
    retval(i) = 0;
  end
end
