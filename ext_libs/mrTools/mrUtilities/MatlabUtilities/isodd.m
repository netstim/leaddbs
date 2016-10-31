% isodd.m
%
%      usage: isodd.m()
%         by: justin gardner
%       date: 10/19/04
%             is number odd
%
function retval = isodd(num)

if (nargin ~= 1)
  help isodd
  return
end


retval =  ~(floor(num/2)*2 == num);

