% stimValsFromNames.m
%
%        $Id:$ 
%      usage: stimVals = stimValsFromNames(stimNames)
%         by: justin gardner
%       date: 10/02/13
%    purpose: Hack to get stimulus values from stimulus names returned form getStimvol
%
function stimVals = stimValsFromNames(stimNames)

% check arguments
if ~any(nargin == [1])
  help stimValsFromNames
  return
end

% go through each stimName
for iVal = 1:length(stimNames)
  thisStimName = stimNames{iVal};
  while ~isempty(thisStimName)
    % parse through looking for space separated stimulus names
    [thisVal thisStimName] = strtok(thisStimName);
    % this is an name=val so break apart by =
    [thisName thisVal] = strtok(thisVal,'=');
    % and get the value
    thisVal = strtok(thisVal,'=');
    % set in output structure
    stimVals(iVal).(thisName) = eval(thisVal);
  end
end

