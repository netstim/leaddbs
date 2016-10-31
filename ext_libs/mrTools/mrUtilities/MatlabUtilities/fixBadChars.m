% fixBadChars.m
%
%      usage: str = fixBadChars(str,<fixList>,<replaceList>,<clipLength>)
%         by: justin gardner
%       date: 04/20/07
%    purpose: takes a string and replaces bad characters not
%    allowed in variable names like space or * with variable name acceptable characters
%             you can also provide your own fixlist, i.e. pairs of
%             things that are the match and replacement, e.g.
%             fixBadChars('remove *this* and replace with that',{'*this*','that'})
% 
%             if you want to mostly keep the fixList, but change one or two item:
%             fixBadChars('change.dot but keep rest',[],{'.','_'});
%  
%             The default is to clip the variable length at 63 characters for a matlab
%             variable, but you can set a different clip length with the fourth argument.
%             Set to 0 for no clipping, [] for default cliping
function str = fixBadChars(str,fixList,replaceList,clipLength)

% check arguments
if ~any(nargin == [1 2 3 4])
  help fixBadChars
  return
end

% this is the list of what characters will map to what
if ieNotDefined('fixList')
  fixList = {{'-','_'},{' ','_'},{'*','star'},{'+','plus'},{'%','percent'},{'[',''},{']',''},{'(',''},{')',''},{'/','_div_'},{'=','_eq_'},{'^','_pow_'},{'.','_period_'},{':','_'},{'&','_and_'},{'!','_bang_'},{'#','_hash_'},{'$','_dollar_'},{'{',''},{'}',''},{'|','_bar_'},{'\','_backslash_'},{';','_'},{'?','_question_'},{',','_comma_'},{'<','_less_'},{'>','_greater_'},{'~','_tilde_'},{'`','_backtick_'}};
  userDefinedFixList = 0;
else
  fixList = cellArray(fixList,2);
  userDefinedFixList = 1;
end

% if we want to keep the default list, but just replace a few characters,
% then the third argument will be a replaceList.
if ~ieNotDefined('replaceList')
  replaceList = cellArray(replaceList,2);
  for i = 1:length(replaceList)
    for j = 1:length(fixList)
      % if we find the one we want to replace
      if isequal(replaceList{i}{1},fixList{j}{1})
	% then replace it
	fixList{j}{2} = replaceList{i}{2};
	break
      end
    end
    % not found, add to our replace list
    fixList{end+1} = replaceList{i};
  end
end

% now swap any occurrences of these characters
for i = 1:length(fixList)
  % look for where we have a bad character
  swaplocs = strfind(str,fixList{i}{1});
  % if any found replace them
  if ~isempty(swaplocs)
    newstr = '';
    swaplocs = [-length(fixList{i}{1})+1 swaplocs];
    for j = 2:length(swaplocs)
      newstr = sprintf('%s%s%s',newstr,str((swaplocs(j-1)+length(fixList{i}{1})):swaplocs(j)-1),fixList{i}{2});
    end
    str = sprintf('%s%s',newstr,str((swaplocs(end)+length(fixList{i}{1})):end));
  end
end

% check for non character at beginning
if ~userDefinedFixList
  if ~isempty(str) && ~isempty(regexp(str,'^[^a-zA-Z]'))
    str = sprintf('x%s',str);
  end
end

if ieNotDefined('clipLength'),clipLength = 63;end
if clipLength > 0
  if length(str) > clipLength
    str = str(1:clipLength);
  end
end

    
