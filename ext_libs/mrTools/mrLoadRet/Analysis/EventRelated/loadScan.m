% loadScan.m
%
%        $Id$
%      usage: loadScan(view,<scanNum>,<groupNum>,<sliceNum>,<precision>,<x>,<y>)
%         by: justin gardner
%       date: 03/20/07
%    purpose: loads a scan into a "d" structure
%             sliceNum is either [] for all slices,
%             slice numbers for slices you want or 0 to not load
%             data
%             precision defaults to 'double' can also be 'single' 
%             x and y are either [] for all rows and cols, or
%             a [min max] for which rows,cols you want
%
function d = loadScan(view,scanNum,groupNum,sliceNum,precision,x,y)

if ~any(nargin == [1 2 3 4 5 6 7])
  help('loadScan');
  return
end
d = [];
if ~isview(view)
  disp(sprintf('(loadScan) First argument is not a view'));
  return
end

% default to loading all slices
if ~exist('sliceNum','var'),sliceNum = [];end
if ~exist('x','var'),x = [];end
if ~exist('y','var'),y = [];end
if ieNotDefined('precision'),precision = 'double';end
if ieNotDefined('scanNum'),scanNum = viewGet(view,'curScan');end
if ieNotDefined('groupNum'),groupNum = viewGet(view,'curGroup');end

% % % % set the group number                      %JB: do not set the group in the view because takes too much time
% % % if ~isempty(groupNum)                       %    and it's not need if no data loading is requested.
% % %   view = viewSet(view,'curGroup',groupNum); %    viewGet can take care of getting info from the right group 
% % % end
% % % groupNum = viewGet(view,'curGroup');
if isempty(groupNum)
  groupNum = viewGet(view,'curGroup');
end


% load parameters
d.ver = 4.5;
d.scanNum = scanNum;
d.groupNum = groupNum;
d.description = viewGet(view,'description',scanNum,groupNum);
d.tr = viewGet(view,'framePeriod',scanNum,groupNum);
d.voxelSize = viewGet(view,'scanvoxelsize',scanNum,groupNum);
d.xform = viewGet(view,'scanSform',scanNum,groupNum);
d.nFrames = viewGet(view,'nFrames',scanNum,groupNum);
d.dim = viewGet(view,'scanDims',scanNum,groupNum);
d.dim(4) = d.nFrames;
d.filename = viewGet(view,'tseriesfile',scanNum,groupNum);
d.filepath = viewGet(view,'tseriespathstr',scanNum,groupNum);
d.expname = getLastDir(fileparts(fileparts(fileparts(d.filepath))));
d.fullpath = fileparts(fileparts(fileparts(fileparts(d.filepath))));
d.volTrigRatio = viewGet(view,'auxParam','volTrigRatio',scanNum,groupNum);

% print out tr for 3d scans to make sure it is right
if viewGet(view,'3D',scanNum,groupNum)
  disp(sprintf('(loadScan) 3D sequence. TR is %0.2f',d.tr));
end

% dispay string to say what we are loading
mrDisp(sprintf('(loadScan) Loading scan %i from group: %s',scanNum,viewGet(view,'groupName',groupNum)));
if length(sliceNum) == 2
  if sliceNum(1) ~= sliceNum(2)
    mrDisp(sprintf(' slices=%i:%i of %i',sliceNum(1),sliceNum(2),viewGet(view,'nSlices',scanNum,groupNum)));
  else
    mrDisp(sprintf(' slice=%i of %i',sliceNum(1),viewGet(view,'nSlices',scanNum,groupNum)));
  end
elseif length(sliceNum) == 1
  mrDisp(sprintf(' slice=%i of %i',sliceNum(1),viewGet(view,'nSlices',scanNum,groupNum)));
end
if length(x) == 2
  if (x(1) ~= 1) || (x(2) ~= d.dim(1))
    mrDisp(sprintf(' x=%i:%i of %i',x(1),x(2),d.dim(1)));
  elseif (x(1) == x(2))
    mrDisp(sprintf(' x=%i of %i',x(1),d.dim(1)));
  end      
elseif length(x) == 1
  mrDisp(sprintf(' x=%i of %i',x(1),d.dim(1)));
end
if length(y) == 2
  if (y(1) ~= 1) || (y(2) ~= d.dim(2))
    mrDisp(sprintf(' y=%i:%i of %i',y(1),y(2),d.dim(2)));
  end
elseif length(y) == 1
  mrDisp(sprintf(' y=%i of %i',y(1),d.dim(2)));
end
mrDisp(sprintf('\n'));

% load the data
if ~isequal(sliceNum,0)
  d.data = loadTSeries(view,scanNum,sliceNum,[],x,y,precision);
else
  d.data = [];
end
	
% set the size of data appropriately
if ~isempty(d.data)
  d.dim = size(d.data);
end

% Dump junk frames
junkFrames = viewGet(view,'junkframes',scanNum,groupNum);
if ~isempty(d.data)
  d.data = d.data(:,:,:,junkFrames+1:junkFrames+d.nFrames);
  % adjust number of volumes to account for the frames that
  % we have just junked
  d.dim(4) = size(d.data,4);
end

% junk frames total is used by getStimvol to adjust
% volumes according to how many volumes have been thrown out
d.junkFrames = viewGet(view,'totalJunkedFrames',scanNum,groupNum);
% we need to add the junk frames junked here to the first
% of the array 'totalJunkedFrames' (one for each scan), or
% if that array is empty then we only have to condiser the
% first junkFrames
if isempty(d.junkFrames)
  d.junkFrames = junkFrames;
else
  d.junkFrames(1) = d.junkFrames(1)+junkFrames;
end
% load dicom header
d.dicom = viewGet(view,'dicom',scanNum,groupNum);

% load stimfile and set traces
d.stimfile = viewGet(view,'stimfile',scanNum,groupNum);

if length(d.junkFrames) ~= length(d.stimfile)
  if (d.junkFrames == 0)
    disp(sprintf('(loadScan) Setting total junk frames for all scans to 0'));
    d.junkFrames = zeros(1,length(d.stimfile));
  end
end

% get any concat info
d.concatInfo = viewGet(view,'concatInfo',scanNum,groupNum);

