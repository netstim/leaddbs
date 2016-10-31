% mlrIsRunning
%
%        $Id$ 
%      usage: tf = mlrIsRunning(verbose)
%         by: justin gardner
%       date: 03/05/10
%    purpose: Retruns whether MLR is running or not
%
function isRunning = mlrIsRunning(verbose)

% check arguments
if ~any(nargin == [0 1])
  help mlrIsRunning
  return
end

if nargin == 0, verbose = 1;end
% init variables
isRunning = 0;
numViews = 0;
mrGlobals;

% go search for views in global
if isfield(MLR,'views')
  for i = 1:length(MLR.views)
    if ~isempty(MLR.views{i}) 
      homeDir = viewGet(MLR.views{i},'homeDir');
      numViews = numViews + 1;
      isRunning = 1;
    end
  end
end

if verbose
  if isRunning
    disp(sprintf('(mlrIsRunning) %i views are open to %s',numViews,homeDir));
  else
    disp(sprintf('(mlrIsRunning) No open views'));
  end
end


