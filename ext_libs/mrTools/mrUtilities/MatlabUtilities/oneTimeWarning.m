% oneTimeWarning.m
%
%        $Id$
%      usage: h = oneTimeWarning(warningName,warningText,<justPrint>)
%         by: justin gardner
%       date: 02/26/08
%    purpose: Brought over from Shani's viewGet function. This
%             will only warn once about something, you give it
%             a warningName and a warningText e.g.:
%
%             oneTimeWarning('someWarning','Something has happened');
%
%             The first time this happens it prints out the
%             warning, after that it won't do anything.
%
%             If the warning is actually just a message, and
%             not a warning (i.e. you just want it printed
%             as is to the matlab buffer), set the justPrint
%             argument to 1.
%
%             If you want the first warning to be an mrWarnDlg
%             and subsequent warnings be printed to the matlab buffer
%             set justPrint = -1;
%
%             To reset the warning so that it will dislay again, do
%             oneTimeWarning('someWarning',0);
%
%             If you want to make this modal, then do the following
%             h = oneTimeWarning('someWarning','I am warning you');
%             if ~isempty(h),uiwait(h);end
%
function h = oneTimeWarning(fieldCheck,warnText,justPrint)

h = [];
% check arguments
if ~any(nargin == [2 3])
  help oneTimeWarning
  return
end

if ieNotDefined('justPrint')
  % default to justPrint if verbose setting is yes
  if isequal(mrGetPref('verbose'),'Yes')
    justPrint = false;
  else
    justPrint = true;
  end
end

% get the warning variable
global gMLRWarning

% make sure the field check name is a valid field name
fieldCheck = fixBadChars(fieldCheck);

% reset warning, if called with a number
if ~isstr(warnText)
  gMLRWarning.(fieldCheck) = [];
  return
end

% if the warning field is not set then...
if ~isfield(gMLRWarning,fieldCheck) | isempty(gMLRWarning.(fieldCheck))
  % set the field to 1 so that it won't print out the next time
  gMLRWarning.(fieldCheck) = 1;
  % and print out a warning
  if justPrint ~= 1
    h = mrWarnDlg(warnText,'Yes');
  else
    disp(warnText);
  end
elseif justPrint == -1
  disp(warnText);
end





