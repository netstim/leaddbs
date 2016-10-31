% iseven.m
%
%      usage: iseven.m()
%         by: justin gardner
%       date: 10/19/04
%             is number even
%
function retval = iseven(num)

if (nargin ~= 1)
  help iseven;
  return
end

retval = floor(num/2)*2 == num;

