% stimconv.m
%
%      usage: stimconv(sig,n)
%         by: justin gardner
%       date: 06/03/05
%
function cmatrix = stimconv(sig,n)

if (nargin ~= 2)
  help stimconv;
  return
end

len = length(sig);
for i=1:n
  cmatrix(1:len,i) = [zeros(1,i-1) sig(1:len-i+1)]';
end
