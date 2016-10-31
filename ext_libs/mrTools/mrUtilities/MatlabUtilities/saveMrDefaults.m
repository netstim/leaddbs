% saveMrDefaults
%
%      usage: saveMrDefaults()
%         by: justin gardner
%       date: 03/17/07
%    purpose: save default positions
%
% Note that the mrDefaults file is usually saved in ~/.mrDefaults
% but that location can be overridden (see mrDefaultsFilename.m)
%
function retval = saveMrDefaults()

% check arguments
if ~any(nargin == [0])
    help saveMrDefaults
    return
end

% get globals
global mrDEFAULTS;

% save figloc
if isfield(mrDEFAULTS,'figloc')
    figloc = mrDEFAULTS.figloc;
else
    figloc = [];
end

if isfield(mrDEFAULTS,'prefs')
    prefs = mrDEFAULTS.prefs;
else
    prefs = [];
end

% check for strange looking prefs before saving
if isempty(prefs) | (length(fields(prefs)) <= 5)
  mrWarnDlg('(saveMrDefaults) Preference variable mrDEFAULTS.prefs does not appear to have any preferences set. Save cancelled');
  return
end
eval(sprintf('save %s figloc prefs -V6;',mrDefaultsFilename));

