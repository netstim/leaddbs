% callbackEval.m
%
%      usage: <output> = callbackEval(callback,<firstArgs>,<lastArgs>)
%        $Id$
%         by: julien besle 
%       date: 16/12/2010
%    purpose: like feval but also handles callbacks in the cell array form
%             firstArgs and lastArgs must be embedded in a cell array 
%             and are added as arguments before and after the arguments 
%             from the callback cell array 
%

function varargout = callbackEval(callback,firstArgs,lastArgs)

if ~ieNotDefined('firstArgs')
  if ~iscell(firstArgs)
    firstArgs = {firstArgs};
  end
  arguments = firstArgs;
else
  arguments = [];
end

if iscell(callback)
  passedArgs = callback(2:end);
  callback = callback{1};
  for iArg = 1:length(passedArgs)
    arguments{end+1} = passedArgs{iArg};
  end 
end

if ~ieNotDefined('lastArgs')
  if ~iscell(firstArgs)
    lastArgs = {lastArgs};
  end
  arguments = [arguments lastArgs];
end

% create the string to call the function
if nargout>=1
  funcall = '[';
  for i=1:nargout
    funcall = sprintf('%svarargout{%i},',funcall,i);
  end
  funcall(end:end+1)= ']=';
else
  funcall = '';
end
funcall = [funcall 'feval(callback'];
for i = 1:length(arguments)
  funcall = sprintf('%s,arguments{%i}',funcall,i);
end
funcall = sprintf('%s);',funcall);

% and call it
eval(funcall);