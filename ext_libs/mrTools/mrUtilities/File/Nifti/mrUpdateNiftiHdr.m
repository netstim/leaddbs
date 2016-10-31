% mrUpdateNiftiHdr.m
%
%      usage: mrUpdateNiftiHdr(<useIdentity>)
%         by: eli merriam
%       date: 03/20/07
%    purpose: Update all nifti headers in groups with
%             the nifti header from the files.
%
%             mrUpdateNiftiHdr;
% 
%             if useIdentity is set, sets all qform and sform
%             matrices to the identity matrix (does not change
%             the actual qform and sform in the nifti hdr, just
%             in the groups, so you can always go back)
%
%             mrUpdateNiftihdr(1);
%
function retval = mrUpdateNiftiHdr(useIdentity)

% check arguments
if ~any(nargin == [0 1])
  help mrUpdateNiftiHdr
  return
end

% default is to update nifti hdrs not set them to identity
if ~exist('useIdentity','var'),useIdentity = 0;,end

% set variables
view = newView;
mrGlobals
updateHdr = 0;

% go through this loop twice, first time is just to tell
% user what is going to change
for passNum = 1:2
  for iGroup = 1:viewGet(view, 'numberofGroups')
    for iScan = 1:viewGet(view, 'nScans', iGroup);
      % load the nifti header from the mrSession file
      curhdr = viewGet(view, 'niftiHdr', iScan, iGroup);

      % if actually going to update nifti headers
      if ~useIdentity
	% load the acutal Nifti header
	filename = viewGet(view, 'tseriesPath', iScan, iGroup);
	hdr = cbiReadNiftiHeader(filename);

	% compare them and update if mismatch
	if ~isequal(curhdr,hdr)
	  updateHdr = 1;
	  % first pass, just display string if something has changed
	  if passNum==1
	    if max(abs(curhdr.sform44(:)-hdr.sform44(:))) > 0.00001
	      disp(sprintf('Nifti hdr for scan %i in %s group has been modified', iScan, viewGet(view, 'groupName',iGroup)));
	    else
	      disp(sprintf('Nifti hdr for scan %i in %s group has been modified (but difference is less than 0.00001)', iScan, viewGet(view, 'groupName',iGroup)));
	    end
	    % second pass, actually do it.
	  else
	    disp(sprintf('Updating nifti hdr for scan %i in %s group', iScan, viewGet(view, 'groupName',iGroup)));
	    view = viewSet(view,'niftiHdr',hdr,iScan,iGroup);
	  end
	end
      % if setting to identity
      else
	updateHdr = 1;
        % first pass, just display string if something has changed
	if passNum == 1
	  disp(sprintf('Nifti hdr qform and sform for scan %i in %s group will be set to identity', iScan, viewGet(view, 'groupName',iGroup)));
	  % second pass, actually do it.
	else
	  disp(sprintf('Setting nifti hdr qform and sform for scan %i in %s group to identity', iScan, viewGet(view, 'groupName',iGroup)));
	  curhdr.qform_code = 1;
	  curhdr.sform_code = 1;
	  curhdr.qform44 = eye(4);
	  curhdr.sform44 = eye(4);
	  view = viewSet(view,'niftiHdr',curhdr,iScan,iGroup);
	end
      end
    end
  end
  % if there is nothing to do, or the user refuses, then return 
  if ~updateHdr
    disp('Headers appear up-to-date');
    return
  elseif ((passNum==1) && ~askuser('Ok to change nifti headers in mrSession?'))
    return
  end
end


% save the mrSession file
saveSession;

% also remove any base anatomies from mrLastView if it is
% there since those might have a different sform
if isfile('mrLastView.mat')
  disp(sprintf('(mrUpdateNiftiHdr) Removing base anatomies from mrLastView'));
  [view viewSettings] = mlrLoadLastView;
  view.baseVolumes = [];
  save mrLastView view viewSettings
end


