% rmParamsField.m
%
%      usage: params = rmParamsField(params)
%         by: julien besle 
%       date: 04/12/2010
%    purpose: rm field from parameters structure output by mrParamsDefault/Dialog
%              $Id: getStimvol.m 1891 2010-11-25 08:06:27Z julien $

function params = mrParamsRemoveField(params,field)

params = rmfield(params,field);
for iField = 1:length(params.paramInfo)
  if strcmp(params.paramInfo{iField}{1},field)
    params.paramInfo(iField) = [];
    break;
  end
end
  
