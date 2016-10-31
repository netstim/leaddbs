% saveSform.m
%
%        $Id$
%      usage: saveSform(sform,[basesformcode],[vol2tal],[vol2mag])
%         by: justin gardner
%       date: 06/20/07
%    purpose: saves sform to mrLoadRet-4.5 mrSession variable and
%             changes header in epis
%
function retval = saveSform(sform,basesformcode,vol2tal,vol2mag)

% check arguments
if nargin == 1
  basesformcode = 1;
  vol2tal = [];
  vol2mag = [];
elseif nargin == 2
  vol2tal = [];
  vol2mag = [];
elseif ~(nargin == 4)
  help saveSform
  return
end

% first check this directory for mrSession.m
path = '';
if isfile('mrSession.mat')
  % if there is a session then ask if the user wants to export to this directory
  answer = questdlg(sprintf('Save alignment to %s?',getLastDir(pwd)),'Export');
  if strcmp(answer,'Cancel'),return,end
  if strcmp(answer,'Yes')
    path = pwd;
  end
end

% if path is not set then ask the user to choose a path
if isempty(path)
  [filename path] = uigetfile('*.mat','Choose a mrSession.mat file to save alignment to');
  if filename == 0,return,end
  if ~strcmp(filename,'mrSession.mat')
    mrWarnDlg(sprintf('(saveSform) %s is not an mrSession file',filename));
    return
  end
end

% start a view in the corresponding location
cd(path);
v = newView;
mrGlobals
% just make sure that the home dir matches
if ~strcmp(MLR.homeDir,path)
  answer = questdlg(sprintf('mrLoadRet is open on session %s? Ok to close and open on %s',getLastDir(MLR.homeDir),getLastDir(path)));
  if ~strcmp(answer,'Yes')
    mrWarnDlg(sprintf('(saveSform) Could not open a view to %s',getLastDir(path)));
    deleteView(v);
    return
  end
  % clear MLR and start over
  deleteView(v);
  clear global MLR;
  v = newView;
end

% ask the user which groups to export to
groupNames = viewGet(v,'groupNames');
for i = 1:length(groupNames)
  paramsInfo{i}{1} = groupNames{i};
  paramsInfo{i}{2} = 1;
  paramsInfo{i}{3} = 'type=checkbox';
  paramsInfo{i}{4} = sprintf('Export alignment to group %s',groupNames{i});
end
params = mrParamsDialog(paramsInfo,'Choose groups to export alignment to');
if isempty(params),return,end

% now go through and update the headers
for iGroup = 1:viewGet(v, 'numberofGroups')
  % see if we are supposed to update the group
  if params.(fixBadChars(viewGet(v,'groupName',iGroup))) %JB: Here we have to fix the bad characters of the group name the same way mrParamsDialog does it
    for iScan = 1:viewGet(v, 'nScans', iGroup);
      % load the nifti header from the mrSession file
      curhdr = viewGet(v, 'niftiHdr', iScan, iGroup);

      % check to see if it is a change (this is just to print user info)
      if ~isfield(curhdr,'sform44') || ~isequal(curhdr.sform44,sform)
	disp(sprintf('Nifti hdr for scan %i in %s group has been modified', iScan, viewGet(v, 'groupName',iGroup)));
      end

      % get the nifti filename
      filename = viewGet(v, 'tseriesPath', iScan, iGroup);
      % check if it is there
      if isfile(filename)
	% load the header
	hdr = mlrImageReadNiftiHeader(filename);
	% set the sform
	hdr = cbiSetNiftiSform(hdr,sform);
	% if the alignment was to a Tal base, need to set the sform_code correctly
	hdr.sform_code = basesformcode;
	
	% and write it back
	hdr = cbiWriteNiftiHeader(hdr,filename);
	% read it back, (I think there is a slight numerical
	% difference in the sform44 from when it is
	% written to when it is read). This doesn't
	% affect anything, but gives better consistency
	% checking for mrUpdateNiftiHeader
	hdr = mlrImageReadNiftiHeader(filename);
	% now save it in the session
	v = viewSet(v,'niftiHdr',hdr,iScan,iGroup);
        
        % save the vol2tal and vol2mag fields
        % only save them if they're not empty, otherwise leave it unchanged
        if ~isempty(vol2tal)
          % check if there's already a vol2tal there that is being overwritten
          oldVol2tal = viewGet(v,'scanvol2tal',iScan,iGroup);
          if ~isempty(oldVol2tal)
            % since aligning to the volume, need to overwrite with the new vol2tal, but warn:
            mrWarnDlg('(saveSform) Overwriting an old Tal xForm');
          end
          v = viewSet(v,'scanvol2tal',vol2tal,iScan,iGroup);
        end
        if ~isempty(vol2mag)
          v = viewSet(v,'scanvol2mag',vol2mag,iScan,iGroup);
        end
      else
	disp(sprintf('(saveSform) Could not open file %s',filename));
      end
    end
  end
end

% save the session
saveSession

% also remove any base anatomies from mrLastView if it is
% there since those might have a different sform
if isfile('mrLastView.mat')
  try
    [view viewSettings] = mlrLoadLastView;
    if isempty(view),return,end
    %here do something something stupid with variable view, otherwise the code analyser is not happy
    %because view is also a function name
    view.baseVolumes;
    %JB: Any base anatomy that is in canonical anatomical base 
    % shouldn't have to be changed, so no need to remove all of them.
    % better to let the user choose wich bases to remove
    %view.baseVolumes = [];
    baselist = selectInList(view,'base','Choose bases to remove',[]);
    if ~isempty(baselist)
      disp(sprintf('(saveSform) Removing base anatomies from mrLastView. Note that if  you have a session open, you still have to remove the bases manually.'));
      view.baseVolumes(baselist) = [];
      if getfield(whos('view'),'bytes')<2e9
	eval(sprintf('save %s view viewSettings -V6;','mrLastView'));
      else
	mrWarnDlg('(mrQuit) Variable view is more than 2Gb, using option -v7.3 to save');
	eval(sprintf('save %s view viewSettings -v7.3;','mrLastView'));
      end
    end
  catch
    disp(sprintf('(saveSform) Error loading mrLastView. This is probably due to some Mathworks issues trying to load a mat file that has old-style figure handles in it. You will need to remove the file mrLastView.mat since Mathworks can no longer seem to load it (complain to them ;-)'));
  end
end
deleteView(v);

