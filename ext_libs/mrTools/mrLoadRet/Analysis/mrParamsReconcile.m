% mrParamsReconcile.m
%
%      usage: mrParamsReconcile(groupName,params,data)
%         by: justin gardner
%       date: 03/13/07
%    purpose: place holder to old file name for old analysis
%             this function should be deleted later (5/21/07)
%             once there are no longer any old analysis that
%             require this function name.
%
function [params data] = mrParamsReconcile(groupName,params,data)

if ieNotDefined('data')
  [params data] = defaultReconcileParams(groupName,params);
else
  [params data] = defaultReconcileParams(groupName,params,data);
end