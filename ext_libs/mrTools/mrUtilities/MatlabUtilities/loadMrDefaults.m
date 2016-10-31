% loadMrDefaults
%
%      usage: loadMrDefaults()
%         by: justin gardner
%       date: 03/17/07
%    purpose: load default positions
%
% Note that the mrDefaults file is usually saved in ~/.mrDefaults
% but that location can be overridden (see mrDefaultsFilename.m)
%
function mrDefaults = loadMrDefaults()

mrDefaults = [];

% check arguments
if ~any(nargin == [0])
  help loadMrDefaults
  return
end

% load the defaults or set them to default values
defaultsFilename = mrDefaultsFilename;
if isfile(defaultsFilename)
  mrDefaults = load(defaultsFilename);
else
    mrDefaults.prefs = [];
    mrDefaults.figloc = [];
end

% set the default values if there are any known fields that have not been set
[prefNames prefDefaults] = mrGetPref;
for i = 1:length(prefNames)
  if ~isfield(mrDefaults.prefs,prefNames{i}) | isempty(mrDefaults.prefs.(prefNames{i}))
    % if we have a cell array, that means a list of different possibilities, default
    % to the top of the list
    if iscell(prefDefaults{i}) && ~isempty(prefDefaults{i})
      mrDefaults.prefs.(prefNames{i}) = prefDefaults{i}{1};
    % otherwise as long as the default is not empty, use it.
    elseif ~isempty(prefDefaults{i})
      mrDefaults.prefs.(prefNames{i}) = prefDefaults{i};
    end
  end
end
  
% check for any figloc that are strange
if ~isempty(mrDefaults.figloc)
  figlocNames = fieldnames(mrDefaults.figloc);
  for i = 1:length(figlocNames)
    % if it isn't of length four then just remove it 
    if length(mrDefaults.figloc.(figlocNames{i})) ~= 4
      mrDefaults.figloc = rmfield(mrDefaults.figloc,figlocNames{i});
    end
  end
end
