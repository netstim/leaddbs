% mrSetNiftiTransforms
%
%      usage: mrSetNiftiTransforms(scanAndGroup,qform44,sform44)
%         by: justin gardner / modified from mrUpdateNiftiHdr
%       date: 03/20/07
%    purpose: sets the qform/sformto the once passed in
%             scanAndGroup sets the nifti header 
%             that should be udpated i.e. if you want to update
%             scans 1 and 4 from group 3 and 2.
%
%             mrSetNiftiTransforms([1 3;4 2],eye(4),sform);
%
%             to leave one transform unchanged do
%
%             mrSetNiftiTransforms([1 3;4 2],eye(4),[]);
%
%
function retval = mrSetNiftiTransforms(scanAndGroup,qform44,sform44)

% check arguments
if ~any(nargin == 3)
  help mrSetNiftiTransforms
  return
end

% set variables
view = newView;
mrGlobals
updateHdr = 0;

% if we get a single number for scanAndGroup it means
% that we want to do all of a certain group
if isequal(size(scanAndGroup),[1 1])
  forceGroup = scanAndGroup;
  scanAndGroup = [];
  scanAndGroup(:,1) = (1:viewGet(view,'nScans',forceGroup))';
  scanAndGroup(:,2) = forceGroup;
end
  
% go through this loop twice, first time is just to tell
% user what is going to change
for passNum = 1:2
  for iGroup = 1:viewGet(view, 'numberofGroups')
    for iScan = 1:viewGet(view, 'nScans', iGroup);
      if ismember([iScan iGroup],scanAndGroup,'rows')
	if passNum == 1
	  % load the nifti header from the mrSession file
	  curhdr = viewGet(view, 'niftiHdr', iScan, iGroup);
	  oldSform = curhdr.sform44;
	  oldQform = curhdr.qform44;
	  disp(sprintf('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'));
	  disp(sprintf('Scan %i in Group %s',iScan,viewGet(view,'groupName',iGroup)));
	  disp(sprintf('++++++++++++++++++++++++ current qform +++++++++++++++++++++++++'));
	  for rownum = 1:4
	    disp(sprintf('%f\t%f\t%f\t%f',oldQform(rownum,1),oldQform(rownum,2),oldQform(rownum,3),oldQform(rownum,4)));
	  end
	  disp(sprintf('++++++++++++++++++++++++ current sform ++++++++++++++++++++++++++'));
	  for rownum = 1:4
	    disp(sprintf('%f\t%f\t%f\t%f',oldSform(rownum,1),oldSform(rownum,2),oldSform(rownum,3),oldSform(rownum,4)));
	  end
	end
        % on second pass, actualy due it
	if passNum == 2
	  % set the transforms
	  if ~isempty(qform44),curhdr.qform44 = qform44;, end
	  if ~isempty(sform44),curhdr.sform44 = sform44;, end
	  % and write the header out
	  %disp(sprintf('Writing %s',sprintf('%s.hdr',stripext(filename))));
	  hdr = cbiWriteNiftiHeader(curhdr,sprintf('%s.hdr',stripext(viewGet(view,'TSeriesPath',iScan,iGroup))));
	  % and change it in the view
	  view = viewSet(view,'niftiHdr',hdr,iScan,iGroup);
	end
      end
    end
  end

  % if there is nothing to do, or the user refuses, then return 
  if ((passNum==1) && ~askuser('Ok to change nifti headers in mrSession?'))
    return
  end
end


% save the mrSession file
saveSession;


