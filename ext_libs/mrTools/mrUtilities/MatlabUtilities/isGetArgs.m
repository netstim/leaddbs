% isGetArgs.m
%
%        $Id:$ 
%      usage: isGetArgs(args)
%         by: justin gardner
%       date: 02/06/12
%    purpose: pass in a cell array to test whether it is a getArgs type argument list
%             specifically, either arguments should be specified as 'arg=x' or be passed in
%             as pairs. If anything is violated returns false.
%
function retval = isGetArgs(args)

% check arguments
if ~any(nargin == [1])
  help isGetArgs
  return
end

% default to false
retval = false;

argNum = 1;

while 1
  % end of list, return true
  if length(args)<argNum
    retval = true;
    return
  end

  % see if the current argument is a string - it should
  % be to specify a name
  if ~isstr(args{argNum}),return,end

  % see if the current argument has an equal sign in which
  % it contains a full string
  eqNum = length(findstr(args{argNum},'='));
  if eqNum > 1,return,end
  % with equal sign we have gotten one complete argument,
  % so continue to next argument
  if eqNum == 1,argNum = argNum+1;continue,end

  % without equal sign the next argument has to be what
  % the value is equal to, so check that we have one
  if length(args)< (argNum+1),return,end
  
  % if it is a string and there is an equal sign, then it
  % is an argument, so should fail
  argNum = argNum+1;
  eqNum = length(findstr(args{argNum},'='));
  if eqNum >= 1,return,end

  % passed, go on to next argument
  argNum = argNum+1;
end


  
