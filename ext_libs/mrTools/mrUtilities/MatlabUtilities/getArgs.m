% evalargs.m
%
%      usage: getArgs(varargin,<validVars>,<verbose=1>)
%         by: justin gardner
%       date: 07/15/08
%    purpose: modification of evalargs to handle validVars
%             better, have a better argument list and evaluate
%             arguments using evalin.
%
%             passed in varargin, creates variables in
%             the calling function determined by the arguments
%             Thus, allows arguments like:
%
%             fun('var1','var2=3','var3',[3 4 5]);
%             will set var1=1, var2=3 and var3=[3 4 5];
%
%             Note that any variable which may not be set in the passed in varargins
%             will need to be declared first in your program (for example):
%
%             function fun(varargin)
%             var1 = [];var2 = [];var3 = [];
%             getArgs(varargin);
%
%             you can have it print
%             out of a list of what is being set by setting the
%             verbose argument to 1.
%
%             getArgs(varargin,[],'verbose=1');
%
%             validVars is a cell array of names that will
%             be checked against to see if the variables that are
%             being set are valid variable names, e.g.
%
%             with a validVars, it will set defaults, and complain
%             if a variable outside of the list is set.
%             getArgs(varargin,{'test1=1','test2=[]','test3=defaultSetting'});
%
%             if you just want to check for valid arguments, but not actually
%             set defaults, then make the validVars list just a list of variable
%             names with no equal signs:
%             getArgs(varargin,{'test1','test2','test3'});
%
%             If you are passing arguments between two functions, one of which will
%             handle one set of variables and the other a different set of values,
%             you can call as follows:
%             [argNames argValues args] = getArgs(varargin,{'fun1var=0'});
%             Then, args will contain an argument list with all the variables
%             *except* fun1var (which will be set in the current function). So 
%             you can then call again in a different function:
%             getArgs(args,{'fun2var=1'});            
%
%             Some other options:
%
%             'verbose=0': Display verbose info
%             'doAssignment=1': Actually set variables.
%             'suppressUnkownArgMessage=0': Set to true to suppress warning for
%                unknown variables
%             'suppressMultipleSettingMessage=false': set to true if you wish to%                suppress message that variables have been set multiple times
%                in the same argument list. If set multiple times, the last
%                one in the list will take precedence
%
function [argNames argValues invalidArgs] = getArgs(args,validVars,varargin)

% check input arguments
if ~any(nargin == [1 2 3 4])
  help getArgs
  return
end

% get our own arguments
if ~ieNotDefined('varargin')
  getArgs(varargin);
end
if ieNotDefined('verbose'),verbose=false;end
if ieNotDefined('doAssignment'),doAssignment=true;end
if ieNotDefined('suppressUnknownArgMessage'),suppressUnknownArgMessage=false;end
if ieNotDefined('suppressMultipleSettingMessage'),suppressMultipleSettingMessage=false;end
% now deal with validVars list
setValidVars = 0;
if ~ieNotDefined('validVars')
  validVars = cellArray(validVars);
  % check to see if we need to set valid values, we only don't 
  % have to do this if the validVars list is all strings, with no equal signs
  for i = 1:length(validVars)
    if ~isstr(validVars{i}) || ~isempty(strfind(validVars{i},'='))
      setValidVars = 1;
    end
  end
  % split into validVarNames and values pairs
  if setValidVars
    [validVarNames validVarValues] = getArgs(validVars,[],'doAssignment=0');
  else
    validVarNames = validVars;
  end
end

% get function name
st = dbstack;
stackNum = 1;
while((stackNum < length(st)) && strcmp(st(stackNum).name,'getArgs'))
  stackNum = stackNum+1;
end
funname = st(stackNum).name;

% loop through arguments
skipnext = 0;argNames = {};argValues = {};invalidArgs={};
for i = 1:length(args)
  % skip if called for
  if skipnext
    skipnext = 0;
    continue
  end
  % evaluate anything that has an equal sign in it
  if isstr(args{i}) && ~isempty(strfind(args{i},'='))
    % if there is nothing after the = sign, then just set to []
    if length(args{i}) == strfind(args{i},'=')
      argNames{end+1} = args{i}(1:strfind(args{i},'=')-1);
      argValues{end+1} = [];
    % if the argument is a numeric, than just set it
    elseif ~isempty(mrStr2num(args{i}(strfind(args{i},'=')+1:end)))
      argNames{end+1} = args{i}(1:strfind(args{i},'=')-1);
      argValues{end+1} = mrStr2num(args{i}(strfind(args{i},'=')+1:end));
    % same for a quoted string
    elseif args{i}(strfind(args{i},'=')+1)==''''
      argNames{end+1} = args{i}(1:strfind(args{i},'=')-1);
      argValues{end+1} = eval(args{i}(strfind(args{i},'=')+1:end));
    % otherwise, we got an unquoted string
    else      
      argNames{end+1} = args{i}(1:strfind(args{i},'=')-1);
      argValues{end+1} = args{i}(strfind(args{i},'=')+1:end);
      % make sure it is not '[]'
      if strcmp(argValues{end},'[]') || strcmp(argValues{end},'{}')
	argValues{end} = [];
      end
    end
  % if it is not evaluated then set it to the next argument,
  % unless there is no next argument in which case set it to 1
  % or if next argument has an equal sign in it
  elseif isstr(args{i}) 
    if (length(args) >= (i+1)) && (~isstr(args{i+1}) || isempty(strfind(args{i+1},'=')))
      % set the variable to the next argument
      argNames{end+1} = args{i};
      argValues{end+1} = args{i+1};
      skipnext = 1;
    else
      % just set the variable to one, since the next argument
      % does not contain a non string
      argNames{end+1} = args{i};
      argValues{end+1} = 1;
    end
  else
    disp(sprintf('(getArgs:%s) Argument %i is not a variable name',funname,i));
  end
end

% test to make sure each variable got set once
[uniqueArgNames uniqueLocs uniqueOriginalLocs] = unique(argNames);
if length(uniqueArgNames) ~= length(argNames)
  for i = 1:length(uniqueLocs)
    varnameCount = sum(uniqueOriginalLocs==i);
    if varnameCount > 1
      if ~suppressMultipleSettingMessage
	disp(sprintf('(%s) Warning: Variable ''%s'' has been set %i times',funname,uniqueArgNames{i},varnameCount));
      end
    end
  end
end

% now go and assign the values
for i = 1:length(argNames)
  % check against variable list
  if ~ieNotDefined('validVars')
    if ~any(strcmp(argNames{i},validVarNames))
      % see if it is just because the case does not match
      caseInsensitiveMatch = find(strcmp(lower(argNames{i}),lower(validVarNames)));
      % if so, give a warning and use the correct capitalization
      if ~isempty(caseInsensitiveMatch)
	disp(sprintf('(%s) Mis-capitalized argument %s changed to %s',funname,argNames{i},validVarNames{first(caseInsensitiveMatch)}));
	argNames{i} = validVarNames{caseInsensitiveMatch};
      % otherwise, give an unknown argument warning
      else
	% if we are not accepting the invalid arguments, then print warning message
	if nargout < 3
	  if ~suppressUnknownArgMessage
	    disp(sprintf('(%s) Unknown argument %s',funname,argNames{i}));
	  end
	end
	invalidArgs{end+1} = argNames{i};
	invalidArgs{end+1} = argValues{i};
      end		
    end
  end
  % assign the variable in caller
  if doAssignment
    assignin('caller',fixBadChars(argNames{i}),argValues{i});
  end
  % if verbose, display what we are doing
  if verbose,dispVarSetting(funname,argNames{i},argValues{i}),end
end

% if we have decided to set the default arguments, then do that
if setValidVars
  for i = 1:length(validVarNames)
    matchingArgNum = find(strcmp(validVarNames{i},argNames));
    if length(matchingArgNum)>1
      matchingArgNum = matchingArgNum(1);
    end
    % if we haven't set the variable, or the variable has been set to []
    if isempty(matchingArgNum) || isempty(argValues{matchingArgNum})
      % assign the default value
      assignin('caller',validVarNames{i},validVarValues{i});
      % if verbose, display what we are doing
      if verbose,dispVarSetting(funname,validVarNames{i},validVarValues{i}),end
    end
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%
%%   dispVarSetting   %%
%%%%%%%%%%%%%%%%%%%%%%%%
function dispVarSetting(funname,varName,varValue)

if isstr(varValue)
  disp(sprintf('(%s) Setting %s=%s',funname,varName,varValue));
elseif isempty(varValue)
  disp(sprintf('(%s) Setting %s=[]',funname,varName));
elseif isnumeric(varValue)
  disp(sprintf('(%s) Setting %s=%s',funname,varName,num2str(varValue)));
else
  disp(sprintf('(%s) Setting %s',funname,varName));
end
