% mlrDispElapsedTime
%
%      usage: str = mlrDispElapsedTime(t)
%         by: justin gardner
%       date: 02/04/04
%       e.g.: mlrDispElapsedTime(130)
%    purpose: t should be an elapsed time in seconds
%
function formatStr = mlrDispElapsedTime(t,varargin)

if nargin < 1
  help mlrDispElapsedTime
  return
end

% get hours,minutes,seconds, ms
hours = floor(t/(60*60));
minutes = floor((t-hours*60*60)/60);
seconds = floor(t-hours*60*60-minutes*60);
ms = floor((t-hours*60*60-minutes*60-seconds)*1000);


% format the string
formatStr = '';
if hours > 0
  formatStr = sprintf('%i hours ',hours);
end
if minutes > 0
  formatStr = sprintf('%s%i min ',formatStr,minutes);
end
if seconds > 0
  formatStr = sprintf('%s%i sec ',formatStr,seconds);
end
formatStr = sprintf('%s%i ms',formatStr,ms);

