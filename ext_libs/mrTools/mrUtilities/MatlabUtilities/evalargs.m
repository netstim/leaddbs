% evalargs.m
%
%      usage: evalargs(varargin,<alwaysUseNextArg>,<verbose>,<validVarnameList>)
%         by: justin gardner
%       date: 12/13/05
%    purpose: passed in varargin, returns a string
%             that once evaluated sets the variables
%             called for. Thus, allows arguments like:
%
%             fun('var1','var2=3','var3',[3 4 5]);
%             will set var1=1, var2=3 and var3=[3 4 5];
%
%             should be run in the folling way:
%
%             function fun(varargin)
%             eval(evalargs(varargin));
%
%             if alwaysUseNextArg is set, then instead of
%             interpreting
%             fun('var1','str')
%               var1=1 str=1
%             will do
%               var1='str'
%
%             when eval is called, you can have it print
%             out of a list of what is being set by setting the
%             verbose argument to 1.
%
%             validVarnameList is a cell array of names that will
%             be checked against to see if the variables that are
%             being set are valid variable names, e.g.
%
%             eval(evalargs(varargin,[],[],{'test1','test2'}));
%             will print out a warning if anything other than test1 and test2 are set
function evalstr = evalargs(args,alwaysUseNextArg,verbose,validVarnameList)

if ~any(nargin == [1 2 3 4])
  evalstr = '';
  help evalargs;
  return
end

% default to always using next argument
if ieNotDefined('alwaysUseNextArg'),alwaysUseNextArg=1;end
if ieNotDefined('verbose'),verbose=0;end

% get function name
st = dbstack;
funname = st(end).name;

% only check the variable names if we are given a list of names
if ieNotDefined('validVarnameList')
  checkValidVarnameList = 0;
else
  checkValidVarnameList = 1;
end

% start the eval string
evalstr = '';
% check arguments in
skipnext = 0;
for i = 1:length(args)
  % skip if called for
  if skipnext
    skipnext = 0;
    continue
  end
  % evaluate anything that has an equal sign in it
  if isstr(args{i}) && ~isempty(strfind(args{i},'='))
    % if the argument is a numeric, than just set it
    if ((exist(args{i}(strfind(args{i},'=')+1:end)) ~= 2) && ...
	~isempty(mrStr2num(args{i}(strfind(args{i},'=')+1:end))))
      evalstr = sprintf('%s%s;',evalstr,args{i});
    % same for a quoted string
    elseif args{i}(strfind(args{i},'=')+1)==''''
      evalstr = sprintf('%s%s;',evalstr,args{i});
    % otherwise, we got a unquoted string, so we need to set the quotes
    else      
      evalstr = sprintf('%s%s''%s'';',evalstr,args{i}(1:strfind(args{i},'=')),args{i}(strfind(args{i},'=')+1:end));
    end
    % any quote needs to be two single quotes
    args{i}(strfind(args{i},''''))='"';
    % if verbose display setting
    if verbose
      evalstr = sprintf('%sdisp(sprintf(''setting: %s''));,',evalstr,args{i});
    end
    % check against validVarnameList
    if checkValidVarnameList && ~any(strcmp(args{i}(1:strfind(args{i},'=')-1),validVarnameList))
      disp(sprintf('(evalargs) Variable %s for function %s is not known',args{i}(1:strfind(args{i},'=')-1),funname));
    end
  % if it is not evaluated then either it means to set the variable
  % or to set the variable to the next argument, we determine this
  % by whether the next argument is a string or not. If it is not
  % a string then it means to set the variable to that argument
  elseif isstr(args{i})
    if (length(args) >= (i+1)) && (~isstr(args{i+1}) || alwaysUseNextArg)
      % set the variable to the next argument
      if ~isstr(args{i+1})
	evalstr = sprintf('%s%s=varargin{%i};',evalstr,args{i},i+1);
        % if verbose display setting
	if verbose
	  evalstr = sprintf('%sdisp(sprintf(''setting: %s=varargin{%i}''));,',evalstr,args{i},i+1);
	end
	if checkValidVarnameList && ~any(strcmp(args{i},validVarnameList))
	  disp(sprintf('(evalargs) Variable %s for function %s is not known',args{i},funname));
	end
      else
	evalstr = sprintf('%s%s=''%s'';',evalstr,args{i},args{i+1});
        % if verbose display setting
	if verbose
	  evalstr = sprintf('%sdisp(sprintf(''setting: %s=''''%s''''''));,',evalstr,args{i},args{i+1});
	end
        % check against validVarnameList
	if checkValidVarnameList && ~any(strcmp(args{i},validVarnameList))
	  disp(sprintf('(evalargs) Variable %s for function %s is not known',args{i},funname));
	end
      end
      skipnext = 1;
    else
      % just set the variable to one, since the next argument
      % does not contain a non string
      evalstr = sprintf('%s%s=1;',evalstr,args{i});
      % if verbose display setting
      if verbose
	evalstr = sprintf('%sdisp(sprintf(''setting: %s=1''));,',evalstr,args{i});
      end
      % check against validVarnameList
      if checkValidVarnameList && ~any(strcmp(args{i},validVarnameList))
	disp(sprintf('(evalargs) Variable %s for function %s is not known',args{i},funname));
      end
    end
  else
    % skip anythin we don't know about
    if ~skipnext
    else
      skipnext = 0;
    end
  end
end

evalstr = sprintf('%s',evalstr);


