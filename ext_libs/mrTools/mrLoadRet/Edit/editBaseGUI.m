% editBaseGUI.m
%
%        $Id$
%      usage: v = editBaseGUI(v)
%         by: justin gardner
%       date: 12/19/07
%    purpose: edit base information
%
function v = editBaseGUI(v)

% check arguments
if ~any(nargin == [1])
  help editBaseGUI
  return
end

% Get current base
baseNum = viewGet(v,'curBase');
if isempty(baseNum)
  mrWarnDlg('(editBaseGUI) No loaded base anatomy');
  return
end

% get associated information
baseName = viewGet(v,'baseName',baseNum);
baseRange = viewGet(v,'baseRange',baseNum);
baseClip = viewGet(v,'baseClip',baseNum);
baseGamma = viewGet(v,'baseGamma',baseNum);
baseCoordMap = viewGet(v,'baseCoordMap',baseNum);

paramsInfo = {};
paramsInfo{end+1} = {'baseName',baseName,'The name of the base anatomy','callback',@editBaseCallback,'callbackArg',v};
paramsInfo{end+1} = {'baseRange',baseRange,'The range of values in the base anatomy','callback',@editBaseCallback,'callbackArg',v};
paramsInfo{end+1} = {'baseMin',baseClip(1),'incdec=[-10 10]','The lower clip value of the base image.','callback',@editBaseCallback,'callbackArg',v};
paramsInfo{end+1} = {'baseMax',baseClip(2),'incdec=[-10 10]','The upper clip value of the base image','callback',@editBaseCallback,'callbackArg',v};
paramsInfo{end+1} = {'baseGamma',baseGamma,'incdec=[-0.1 0.1]','minmax=[0 inf]','The gamma to use to display the base image','callback',@editBaseCallback,'callbackArg',v};
paramsInfo{end+1} = {'setDefaults',0,'type=pushbutton','buttonString','Set defaults','callback',@defaultEditBaseCallback,'callbackArg',v,'passParams=1'};

% deal with coordMap fields
if ~isempty(baseCoordMap)
  baseCoordMapFields = fields(baseCoordMap);
  editableBaseCoordMapFields = {};
  for fieldNum = 1:length(baseCoordMapFields)
    % look for fields that are filenames or directories.
    % these are the ones we allow editing on
    if isstr(baseCoordMap.(baseCoordMapFields{fieldNum}))
      if ~isempty(strfind(baseCoordMapFields{fieldNum},'FileName')) || ~isempty(strfind(baseCoordMapFields{fieldNum},'Dir')) || ~isempty(strfind(baseCoordMapFields{fieldNum},'path'))
	editableBaseCoordMapFields{end+1} = baseCoordMapFields{fieldNum};
	paramsInfo{end+1} = {baseCoordMapFields{fieldNum},baseCoordMap.(baseCoordMapFields{fieldNum}),sprintf('Name of %s from which this surface/flat was generated')};
      end
    end
  end
end


% display parameters
params = mrParamsDialog(paramsInfo,'Edit base info');

% user hit cancel, so return parameters back to what they used to be
if isempty(params)
  v = viewSet(v,'baseName',baseName,baseNum);
  v = viewSet(v,'baseRange',baseRange,baseNum);
  v = viewSet(v,'baseMin',baseClip(1),baseNum);
  v = viewSet(v,'baseMax',baseClip(2),baseNum);
  return
end

% for the baseCoordMap fields, we only set them after user hit ok, not on the fly
if ~isempty(baseCoordMap) && ~isempty(params)
  fieldChange = 0;
  for fieldNum = 1:length(editableBaseCoordMapFields)
    if ~strcmp(baseCoordMap.(editableBaseCoordMapFields{fieldNum}),params.(editableBaseCoordMapFields{fieldNum}))
      baseCoordMap.(editableBaseCoordMapFields{fieldNum}) = params.(editableBaseCoordMapFields{fieldNum});
      fieldChange = 1;
    end
  end
  if fieldChange
    mrWarnDlg('(editBaseGUI) Note that changing the names of surfaces in the base structure does not actually change the surface/flat map. If you want to change the surfaces used, you will need to use Load Flat or Load Surface to create a new flat/surface.');
    viewSet(v,'baseCoordMap',baseCoordMap,baseNum);
  end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   editBaseCallback   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%
function editBaseCallback(v,params)

% get current baseNum
baseNum = viewGet(v,'curBase');

% change info
v = viewSet(v,'baseName',params.baseName,baseNum);
v = viewSet(v,'baseRange',params.baseRange,baseNum);
v = viewSet(v,'baseMin',params.baseMin,baseNum);
v = viewSet(v,'baseMax',params.baseMax,baseNum);
v = viewSet(v,'baseGamma',params.baseGamma,baseNum);

% refresh the display
refreshMLRDisplay(viewGet(v,'viewNum'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   default clipping values taken from isbase   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function clip = defaultEditBaseCallback(v,params)

% get current baseNum
baseNum = viewGet(v,'curBase');

% get base image
b = viewGet(v,'base',baseNum);
image = b.data;

% Choose default clipping based on histogram
histThresh = length(image(:))/1000;
[cnt, val] = hist(image(:),100);
goodVals = find(cnt>histThresh);
clipMin = val(min(goodVals));
clipMax = val(max(goodVals));
clip = [clipMin,clipMax];

v = viewSet(v,'baseMin',clipMin,baseNum);
v = viewSet(v,'baseMax',clipMax,baseNum);
v = viewSet(v,'baseRange',[min(image(:)) max(image(:))]);
v = viewSet(v,'baseGamma',1);

% refresh the display
refreshMLRDisplay(viewGet(v,'viewNum'));

% reset the parameters
params.baseMin = clipMin;
params.baseMax = clipMax;
params.baseRange = [min(image(:)) max(image(:))];
params.baseGamma = 1;
mrParamsSet(params);
