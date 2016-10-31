function preferences = mrSetPref(pref,value,verbose)
%
% preferences = mrSetPref(pref,value)
%
% Replaces Matlab's setpref function. Sets a field in the global variable
% mrDEFAULTS.preferences, which is a structure with fields for each
% preference. Returns the preferences structure.
%
% Examples:
%   mrSetPref('verbose','Yes');
%   mrSetPref('verbose','No');
%   mrSetPref('site','NYU');
%   mrSetPref('niftiFileExtension','.img');
%   mrSetPref('niftiFileExtension','.nii');
%   mrSetPref('interpMethod','nearest');
%      Options: 'nearest','linear','spline','cubic'
%   mrSetPref('overwritePolicy','Ask');
%      Options: 'Ask','Merge','Rename','Overwrite'
%
% Note that the mrDefaults file is usually saved in ~/.mrDefaults
% but that location can be overridden (see mrDefaultsFilename.m)
%
% you can reset a value back to its default value by doing
%
%  mrSetPref('interpMethod');
% 
% For preferences that are registered within the mrGetPref function
% this will do some checking of captilization and default values,
% but if you wish to use a preference value that has not been so
% registered you can do so and suppress warning that the preference
% is unknown:
%
% mrSetPref('unknwonPrefName','prefValue',false);
%
% djh, 5/2007

if ~any(nargin == [1 2 3])
  help mrSetPref;
  return
end

if nargin < 3,verbose = true;end
global mrDEFAULTS

% if mrDEFAULTS is empty, then we should try to load it
if isempty(mrDEFAULTS)
  mrDEFAULTS = loadMrDefaults;
end

% get pref names and defaults
[prefNames prefDefaults] = mrGetPref;

% check for a known preference
prefNum = find(strcmp(lower(pref),lower(prefNames)));

% check if value is not set, if so then reset with default
if ~exist('value','var')
  if ~isempty(prefNum) && ~isempty(prefDefaults{prefNum})
    if iscell(prefDefaults{prefNum}) 
      value = prefDefaults{prefNum}{1};
    else
      value = prefDefaults{prefNum};
    end
  else
    mrWarnDlg(sprintf('(mrSetPref) No default value for preference %s',pref));
    value = [];
  end
end


if isempty(prefNum)
  if verbose
    % print message for unknown preference
    disp(sprintf('(mrSetPref) Unknown preference %s',pref));
  end
else
  % this will fix the caps on the prefs name
  pref = prefNames{prefNum};
  % check for a known default
  if ~isempty(prefDefaults{prefNum}) && iscell(prefDefaults{prefNum}) && ~iscell(prefDefaults{prefNum}{1})
    prefDefaultNum = find(strcmp(lower(value),lower(prefDefaults{prefNum})));
    % print message if it is not known
    if isempty(prefDefaultNum)
      mrWarnDlg(sprintf('(mrSetPref) Value %s for preference %s is not known',value,pref));
    else
      % fix caps
      value = prefDefaults{prefNum}{prefDefaultNum};
    end
  end
end

mrDEFAULTS.prefs = setfield(mrDEFAULTS.prefs,pref,value);
preferences = mrDEFAULTS.prefs;
saveMrDefaults;