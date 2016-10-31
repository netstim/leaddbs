% first.m
%
%      usage: first()
%         by: justin gardner
%       date: 12/09/04
%    purpose: returns first element of input array
function retval = first(x)

if (isempty(x))
  retval = [];
elseif iscell(x)
  retval = x{1};
else
  retval = x(1);
end