% mrParamsDefault.m
%
%      usage: mrParamsDefault(paramsinfo)
%         by: justin gardner
%       date: 03/13/07
%    purpose: returns default params given paramInfo
%             see wiki for details
%
function params = mrParamsDefault(paramInfo)

% check arguments
if ~any(nargin == [1])
  help mrParamsDefault
  return
end

% parse the input parameter string
[vars varinfo] = mrParamsParse(paramInfo);
for i = 1:length(vars)
  % if it is a group, then return the cell array of all values
  if isfield(varinfo{i},'group')
    if strcmp(varinfo{i}.type,'numeric') || strcmp(varinfo{i}.type,'checkbox')
      % if it is something that should be numeric, than make
      % it into an array
      for j = 1:length(varinfo{i}.allValues)
	if isstr(varinfo{i}.allValues{j})
	  params.(varinfo{i}.name)(j) = str2num(varinfo{i}.allValues{j});
	else
	  params.(varinfo{i}.name)(j) = varinfo{i}.allValues{j};
	end
      end
    % check for popupmenu
    elseif strcmp(varinfo{i}.type,'popupmenu')
      % take the top of the popup list
      for j = 1:length(varinfo{i}.allValues)
	params.(varinfo{i}.name){j} = varinfo{i}.allValues{j}{1};
      end
      % otherwise, just copy the list of parameters
    else
      params.(varinfo{i}.name) = varinfo{i}.allValues;
    end
    continue;
  end
  % make sure it is not a contingent value that has been shut
  % off, first get value it is contingent on
  if isfield(varinfo{i},'contingent')
    contingentValue = varinfo{varinfo{i}.contingentOn}.value;
    if isstr(contingentValue),contingentValue = str2num(contingentValue);,end
    % if it has been shut down, give the parameter an empty
    % value and continue on
    if isequal(contingentValue,0)
      params.(varinfo{i}.name) = [];
      continue
    end
  end
  % otherwise set the parameter to the default value
  if ~strcmp(varinfo{i}.type,'popupmenu')
    params.(varinfo{i}.name) = varinfo{i}.value;
  else
    params.(varinfo{i}.name) = varinfo{i}.value{1};
  end
  % change to numeric if need be
  if strcmp(varinfo{i}.type,'numeric') || strcmp(varinfo{i}.type,'checkbox')
    if isstr(params.(varinfo{i}.name))
      params.(varinfo{i}.name) = str2num(params.(varinfo{i}.name));
    end
  end
  % check to see if this is an array that could not be displayed
  if isfield(varinfo{i},'tooBigArrayValue')
    % if so then set to the original array value
    params.(varinfo{i}.name) = varinfo{i}.tooBigArrayValue;
  end
end
params.paramInfo = paramInfo;

