% last.m
%
%      usage: last.m()
%         by: justin gardner
%       date: 02/10/05
%    purpose: get last element of array
%
function retval = last(x)

if isempty(x)
  retval = [];
elseif iscell(x)
  retval = x{end};
else
  retval = x(end);
end

