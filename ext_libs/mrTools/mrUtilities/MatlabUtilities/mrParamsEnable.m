% mrParamsEnable.m
%
%        $Id$ 
%      usage: mrParamsEnable(paramname,enable)
%         by: justin gardner
%       date: 07/23/09
%    purpose: enable or disable a a parameter in an mrParamsDialog by name. Enable should be 1 or 0
%
function retval = mrParamsEnable(paramname,enable)

% check arguments
if ~any(nargin == [2])
  help mrParamsEnable
  return
end

if isequal(enable,1)
  enable = 'on';
elseif isequal(enable,0)
  enable = 'off';
end

global gParams;
if isempty(gParams)
  disp('(mrParamsEnable) mrParamsDialog is not running');
  return
end

for i = 1:length(gParams.varinfo)
  if strcmp(gParams.varinfo{i}.name,paramname)
    % check for contingentOn
    thisEnable = enable;
    if strcmp(enable,'on') && isfield(gParams.varinfo{i},'contingentOn')
      % if it is contingent on something that has a value of 0, then
      % keep it off
      params = mrParamsGet('getUnenabled=1','paramNum',gParams.varinfo{i}.contingentOn);
      controlName = gParams.varinfo{gParams.varinfo{i}.contingentOn}.name;
      if isfield(params,controlName) && isequal(params.(controlName),0)
	thisEnable = 'off';
      end
    end
    for j = 1:length(gParams.ui.varentry{i})
      % set enable on main field
      set(gParams.ui.varentry{i}(j),'enable',thisEnable);
    end
    % set enable on inc/dec
    for j = 1:length(gParams.ui.incdec{i})
      if ~isempty(gParams.ui.incdec{i}{j})
	set(gParams.ui.incdec{i}{j},'enable',thisEnable);
      end
    end
    drawnow
    return
  end
end

disp(sprintf('(mrParamsEnable) Could not find parameter %s',paramname));



