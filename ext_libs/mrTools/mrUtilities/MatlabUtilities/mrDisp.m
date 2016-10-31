% mrDisp.m
%
%      usage: mrDisp(str)
%         by: justin gardner
%       date: 06/03/08
%    purpose: prints w/out newline and with a flush
%
function mrDisp(str)

% if this is being called it means the mex file doesn't exist,
% so just print out (this won't flush though--preventing updating
% text print outs like for disppercent)
fprintf(1,str);

