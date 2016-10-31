% setFramePeriod.m
%
%        $Id: setFramePeriod.m,v 1.3 2007/09/22 11:00:35 justin Exp $
%      usage: setFramePeriod(framePeriod)
%         by: justin gardner
%       date: 06/20/07
%    purpose: saves framePeriod to mrLoadRet-4.5 mrSession variable and
%             changes header in epis (i.e. pixdim(5) to what is passed in)
%             note that what gets saved in the header is set in ms, but
%             the argument passed in should be in seconds. e.g. to change
%             the frame period to 1.5 seconds:
%
%             setFramePeriod(1.5);
%
function retval = setFramePeriod(framePeriod)

% check arguments
if ~any(nargin == [1])
  help setFramePeriod
  return
end

% first check this directory for mrSession.m
path = '';
if isfile('mrSession.mat')
  % if there is a session then ask if the user wants to export to this directory
  answer = questdlg(sprintf('Save frame period to %s?',getLastDir(pwd)),'setFramePeriod');
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
    mrWarnDlg(sprintf('(setFramePeriod) %s is not an mrSession file',filename));
    return
  end
end

% start a view in the corresponding location
cd(path);
v = newView('Volume');
mrGlobals
% just make sure that the home dir matches
if ~strcmp(MLR.homeDir,path)
  answer = questdlg(sprintf('mrLoadRet is open on sesion %s? Ok to close and open on %s',getLastDir(MLR.homeDir),getLastDir(path)));
  if ~strcmp(answer,'Yes')
    mrWarnDlg(sprintf('(setFramePeriod) Could not open a view to %s',getLastDir(path)));
    deleteView(v);
    return
  end
  % clear MLR and start over
  deleteView(v);
  clear global MLR;
  v = newView('Volume');
end

% ask the user which groups to export to
groupNames = viewGet(v,'groupNames');
for i = 1:length(groupNames)
  paramsInfo{i}{1} = groupNames{i};
  paramsInfo{i}{2} = 1;
  paramsInfo{i}{3} = 'type=checkbox';
  paramsInfo{i}{4} = sprintf('Save framePeriod of %0.2f to group %s',framePeriod,groupNames{i});
end
params = mrParamsDialog(paramsInfo,'Choose groups to export alignment to');
if isempty(params),return,end

% now go through and update the headers
for iGroup = 1:viewGet(v, 'numberofGroups')
  % see if we are supposed to update the group
  if params.(fixBadChars(viewGet(v,'groupName',iGroup)))
    for iScan = 1:viewGet(v, 'nScans', iGroup);
      % load the nifti header from the mrSession file
      curhdr = viewGet(v, 'niftiHdr', iScan, iGroup);
      % check to see if it is a change (this is just to print user info)
      if isfield(curhdr,'pixdim') && ~isequal(curhdr.pixdim(5),framePeriod*1000)
	disp(sprintf('Frame period of %0.2f for scan %i in %s group has been set to %0.2f', curhdr.pixdim(5)/1000, iScan, viewGet(v, 'groupName',iGroup),framePeriod));
      end
      % get the nifti filename
      filename = viewGet(v, 'tseriesPath', iScan, iGroup);
      % check if it is there
      if isfile(filename)
	% load the header
	hdr = mlrImageReadNiftiHeader(filename);
	newUnits = bitand(hdr.xyzt_units,7)+16;
	if (newUnits ~= hdr.xyzt_units)
	  disp(sprintf('(setFramePeriod) Changing xyzt_units in nifit header from %s (%i) to %s (%i)',dec2bin(hdr.xyzt_units),hdr.xyzt_units,dec2bin(newUnits),newUnits));
	  hdr.xyzt_units = newUnits;
	end
	% set the frameperiod
	hdr.pixdim(5) = framePeriod*1000;
	% and write it back
	hdr = cbiWriteNiftiHeader(hdr,filename);
	% set the scan params
	scanParams = viewGet(v,'scanParams',iScan,iGroup);
	scanParams.framePeriod = framePeriod;
	% now save it in the session
	viewSet(v,'scanParams',scanParams,iScan,iGroup);
	viewSet(v,'niftiHdr',hdr,iScan,iGroup);
      else
	disp(sprintf('(setFramePeriod) Could not open file %s',filename));
      end
    end
  end
end

% save the session
saveSession

deleteView(v);

