% mrStr2num
%
%        $Id$
%      usage: num = mrStr2n(str)
%         by: justin gardner
%       date: 07/17/07
%    purpose: returns number or empty if string is not a number
%             matlab's str2num is *very* annoying since it
%             evaluates strings, so that if your string happens to
%             have, say a function name in it then that function
%             will be called, instead of just returning []
%
function retval = mrStr2num(str)

retval = [];

% check arguments
if ~any(nargin == [1])
  help mrStr2num
  return
end

% check for empy
if isempty(str)
  retrval = [];
  return;
end

% remove from the string any nan/inf for testing
% since these are valid strings to have.
teststr = fixBadChars(lower(str),{{'nan',''},{'inf',''},{'-',''}});

% now go through the string and check
% to see if any space delimited part of
% it is actually a function
while ~isempty(teststr)
  [tok teststr] = strtok(teststr,' ');
  % if the token is a function then return, this is not a number
  if any(exist(tok) == [3 5])
    return;
  end
  % check to see if it is an m-file, make sure to put the .m on
  % otherwise, matlab will pick up regular files (for some reason
  % I had a file called "1" which caused a lot of problems :-)
  tok = sprintf('%s.m',stripext(tok));
  if any(exist(tok) == [2])
    return;
  end
end

% whatever made it here is not a function,
% so run str2num on it
retval = str2num(str);
