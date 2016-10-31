% quickReconcile.m
%
%      usage: quickReconcile(groupName,params)
%         by: justin gardner
%       date: 04/04/07
%    purpose: function used by quickAnalysis routines to satisfy
%             having a reconcile function
%
function [params data] = quickReconcile(groupName,params,data)

if ieNotDefined('params')
  params = [];
end
if ieNotDefined('data')
  data = [];
end

