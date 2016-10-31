% mrFixHdr
%
%      usage: mrFixNiftiHdr(<forceUpdate>,<copyFrom>)
%         by: justin gardner / modified from mrUpdateNiftiHdr
%       date: 03/20/07
%    purpose: "Fixes" nifti headers, by copying them from the original
%             nifti file they were created from. This is useful for undoing
%             the problem that FSL causes when you run motion correction.
%
%             mrFixHdr;
% 
%             force update is a list of all scan and groups
%             that should be udpated i.e. if you want to update
%             scans 1 and 4 from group 3 and 2.
%
%             mrFixNiftiHdr([1 3;4 2]);
%
%             If you want to copy all headers for a certain scan from
%             the orginal, you can do (e.g. for group 2)
%
%             mrFixNiftiHdr(2);
% 
%             if copy from is set, will copy the nifti header from
%             the scan group, set in copyFrom, so to copy from
%             the header for scan1, group1
%
%             mrFixNiftiHdr([1 3;4 2],[1 1]);
%
%
function retval = mrFixNiftiHdr(forceUpdate,copyFrom)

% check arguments
if ~any(nargin == [0 1 2])
  help mrFixNiftiHdr
  return
end

% default is to update nifti hdrs not set them to identity
if ~exist('forceUpdate','var'),forceUpdate = [];,end
if ~exist('copyFrom','var'),copyFrom = [];,end

% set variables
view = newView;
mrGlobals
updateHdr = 0;

% if we get a single number for forceUpdate it means
% that we want to do all of a certain group
if isequal(size(forceUpdate),[1 1])
  forceGroup = forceUpdate;
  forceUpdate = [];
  forceUpdate(:,1) = (1:viewGet(view,'nScans',forceGroup))';
  forceUpdate(:,2) = forceGroup;
end
  
% set the string to display
if isempty(forceUpdate)
  dispReasonString = 'mismatch with qform44: ';
else
  dispReasonString = '';
end

% go through this loop twice, first time is just to tell
% user what is going to change
for passNum = 1:2
  for iGroup = 1:viewGet(view, 'numberofGroups')
    for iScan = 1:viewGet(view, 'nScans', iGroup);
      % load the nifti header from the mrSession file
      curhdr = viewGet(view, 'niftiHdr', iScan, iGroup);

      % load the acutal Nifti header
      filename = viewGet(view, 'tseriesPath', iScan, iGroup);
      hdr = cbiReadNiftiHeader(filename);
      % get dimensions of voxels from headr
      voxdim = hdr.pixdim(2:4)';
      % get dimensions of voxels from qform44
      xdim = sqrt(sum(hdr.qform44(:,1).^2));
      ydim = sqrt(sum(hdr.qform44(:,2).^2));
      sdim = sqrt(sum(hdr.qform44(:,3).^2));
      voxdimFromQform = [xdim ydim sdim];
      % round to nearest 10000
      voxdim = round(voxdim*1000)/1000;
      voxdimFromQform = round(voxdimFromQform*1000)/1000;
      if (isempty(forceUpdate) && ~isequal(voxdim,voxdimFromQform)) || ismember([iScan iGroup],forceUpdate,'rows')
	% get original header
	[os og] = viewGet(view,'originalScanNum',iScan,iGroup);
	% if we want to copy from a specific one
	if ~isempty(copyFrom)
	  os = copyFrom(1);og = copyFrom(2);
	end
	updateHdr = 1;
	% print out info about what is going on
	if ~isempty(os)
	  disp(sprintf('(%s:%i) pixdim [%s] %s copy header from (%s:%i) to (%s:%i)',viewGet(view,'groupName',iGroup),iScan,num2str(voxdim),dispReasonString,viewGet(view,'groupName',og(1)),os(1),viewGet(view,'groupName',iGroup),iScan));
	else
	  disp(sprintf('(%s:%i) pixdim [%s] %s No original header to copy from',viewGet(view,'groupName',iGroup),iScan,num2str(voxdim),dispReasonString));
	  % user will have to specify which group/scan to copy from
	  if passNum == 2
	    og = getnum('Which group do you want to copy from? ',1:viewGet(v,'numberOfGroups'));
	    os = getnum('Which scan do you want to copy from? ',1:viewGet(v,'numberOfScans',og));
	  end
	end
	% on second pass, actualy due it
	if passNum == 2
	  % get the source header
	  srchdr = viewGet(view,'niftiHdr',os(1),og(1));
	  srcVoxelSize = viewGet(view,'scanVoxelSize',os(1),og(1));
	  % and write the header out
	  %disp(sprintf('Writing %s',sprintf('%s.hdr',stripext(filename))));
	  srchdr = cbiWriteNiftiHeader(srchdr,sprintf('%s.hdr',stripext(filename)));
	  % and change it in the view
	  view = viewSet(view,'niftiHdr',srchdr,iScan,iGroup);
	  scanParams = viewGet(view,'scanParams',iScan,iGroup);
	  scanParams.voxelSize = srcVoxelSize;
	  view = viewSet(view,'scanParams',scanParams,iScan,iGroup);
	end
      else
	disp(sprintf('(%s:%i) pixdim [%s] matches with qform44',viewGet(view,'groupName',iGroup),iScan,num2str(voxdim)));
      end
    end
  end

  % if there is nothing to do, or the user refuses, then return 
  if ~updateHdr
    disp('Headers do not need to be fixed');
    return
  elseif ((passNum==1) && ~askuser('Ok to change nifti headers in mrSession?'))
    return
  end
end


% save the mrSession file
saveSession;


