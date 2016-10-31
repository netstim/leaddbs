% groupInfo.m
%
%      usage: groupInfo(groupNum)
%         by: justin gardner
%       date: 11/08/06
%    purpose: print out information about scans in group
%       e.g.: groupInfo(1);
%             groupInfo('Raw');
%
%             to get info on all groups:
%             groupInfo
% 
%             to print info on groups with scanSforms
%             groupInfo('Raw',1)
%
%             an optional return argument returns a structure
%             with the group information
%             info = groupInfo('Raw');
%
function retval = groupInfo(groupNum,verbose,dispDialog)

% check arguments
if ~any(nargin == [0 1 2 3])
  help groupInfo
  return
end

retval = [];
if ieNotDefined('verbose'),verbose = 0;end
if ieNotDefined('dispDialog'),dispDialog = 0;end
view = newView;
if isempty(view),return,end

% check home dir
homeDir = viewGet(view,'homeDir');
if ~isequal(homeDir,pwd)
  disp(sprintf('(groupInfo) Current directory (%s) is not the home directory of the session',getLastDir(pwd)));
  deleteView(view);
  return
end

%%%%%%%%%%%%%%%%%%%%%%%%%
% groupInfo with no arguments
%%%%%%%%%%%%%%%%%%%%%%%%%
% if groupNum was not given then display info about all groups
if ieNotDefined('groupNum')
  % get session info
  subject = viewGet(view,'subject');
  description = viewGet(view,'sessionDescription');
  magnet = viewGet(view,'magnet');
  operator = viewGet(view,'operator');
  coil = viewGet(view,'coil');
  protocol = viewGet(view,'protocol');
  % display session info
  disp(sprintf('%s',repmat('=',1,40)));
  disp(sprintf('homeDir: %s',homeDir));
  disp(sprintf('description: %s',description));
  disp(sprintf('operator: %s subject: %s',operator,subject));
  disp(sprintf('magnet: %s coil: %s protocol: %s',magnet,coil,protocol));
  disp(sprintf('%s',repmat('=',1,40)));
  % display each group
  for g = 1:viewGet(view,'nGroups')
    % get group info
    retval(g).groupName = viewGet(view,'groupName',g);
    retval(g).numScans = viewGet(view,'nScans',g);
    % use du to get disk usage
    [status result] = system(sprintf('du -k -d 0 %s',viewGet(view,'datadir',g)));
    retval(g).dirSize = str2num(strtok(result));
    if (retval(g).dirSize > 1000000)
      dirSize = sprintf('%0.1fG',retval(g).dirSize/1000000);
    elseif (retval(g).dirSize > 1000)
      dirSize = sprintf('%0.1fM',retval(g).dirSize/1000);
    else
      dirSize = sprintf('%iK',retval(g).dirSize);
    end
    % display group info
    disp(sprintf('%i: %s (%i scans) %s',g,retval(g).groupName,retval(g).numScans,dirSize));
  end
  deleteView(view);
  % don't return anything unless asked for
  if nargout == 0,clear retval,end
  return
end

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

% now go through scans and print information
paramsInfo = {};
for s = 1:viewGet(view,'numberOfScans',groupNum)
  % grab info
  retval(s).description = viewGet(view,'description',s,groupNum);
  retval(s).scanVoxelSize = viewGet(view,'scanVoxelSize',s,groupNum);
  retval(s).tr = viewGet(view,'framePeriod',s,groupNum);
  retval(s).totalFrames = viewGet(view,'totalFrames',s,groupNum);
  retval(s).filename = viewGet(view,'tSeriesFile',s,groupNum);
  retval(s).originalFilename = viewGet(view,'originalFilename',s,groupNum);
  retval(s).originalGroupname = viewGet(view,'originalGroupname',s,groupNum);
  retval(s).stimFilename = viewGet(view,'stimFilename',s,groupNum);
  retval(s).scanDims = viewGet(view,'scanDims',s,groupNum);
  retval(s).totalJunkedFrames = viewGet(view,'totalJunkedFrames',s,groupNum);
  retval(s).junkFrames = viewGet(view,'junkFrames',s,groupNum);
  
  paramsInfo{end+1} = {sprintf('Scan%02i',s),sprintf('%s: %i volumes',retval(s).description,retval(s).totalFrames),'editable=0'};
  
  % display info
  disp(sprintf('%i: %s',s,retval(s).description));
  disp(sprintf('   Filename: %s GroupName: %s',retval(s).filename,groupName));
  for i = 1:length(retval(s).originalFilename)
    disp(sprintf('   Original Filename: %s Group: %s',retval(s).originalFilename{i},retval(s).originalGroupname{i}));
  end
  for i = 1:length(retval(s).stimFilename)
    if iscell(retval(s).stimFilename{i})
      for j = 1:length(retval(s).stimFilename{i})
	disp(sprintf('   StimFilename: %s',retval(s).stimFilename{i}{j}));
      end
    else
      disp(sprintf('   StimFilename: %s',retval(s).stimFilename{i}));
    end
  end

  disp(sprintf('   junkFrames=[%s] totalJunkedFrames=[%s]',num2str(retval(s).junkFrames),num2str(retval(s).totalJunkedFrames)));
  disp(sprintf('   voxelSize=[%0.1f %0.1f %0.1f] TR=%0.4f Dims: [%i %i %i] Volumes=%i',retval(s).scanVoxelSize(1),retval(s).scanVoxelSize(2),retval(s).scanVoxelSize(3),retval(s).tr,retval(s).scanDims(1),retval(s).scanDims(2),retval(s).scanDims(3),retval(s).totalFrames));

  % if verbose is set to 1 or above, then show scan transform
  if verbose >= 1 
    retval(s).scanSform = viewGet(view,'scanSform',s,groupNum);
    if ~isempty(retval(s).scanSform)
      disp(sprintf('   scanSform:'));
      for i = 1:4
	disp(sprintf('   %0.3f %0.3f %0.3f %0.3f',retval(s).scanSform(i,1),retval(s).scanSform(i,2),retval(s).scanSform(i,3),retval(s).scanSform(i,4)));
      end
    else
      disp(sprintf('   scanSform: []'));
    end
  end
end

if dispDialog
  mrParamsDialog(paramsInfo,'groupInfo',1.5);
end
deleteView(view);

% don't return anything unless asked for
if nargout == 0,clear retval,end
