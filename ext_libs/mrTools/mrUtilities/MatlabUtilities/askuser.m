% askuser.m
%
%      usage: askuser(question,<toall>,<useDialog>)
%         by: justin gardner
%       date: 02/08/06
%    purpose: ask the user a yes/no question. Question is a string or a cell arry of strings with the question. This
%             function will return 1 for yes and 0 for no. If toall is set to 1, then
%             'Yes to all' will be an option, which if selected will return inf
%
function r = askuser(question,toall,useDialog)

% check arguments
if ~any(nargin == [1 2 3])
  help askuser
  return
end

if ieNotDefined('toall'),toall = 0;,end
if ieNotDefined('useDialog'),useDialog=0;end

% if explicitly set to useDialog then use dialog, otherwise
% check verbose setting
if useDialog
  verbose = 1;
else
  verbose = mrGetPref('verbose');
  if strcmp(verbose,'Yes'),verbose = 1;else,verbose = 0;end
end

r = [];

question=cellstr(question); %convert question into a cell array
  

while isempty(r)
  % ask the question
  if ~verbose
    % not verbose, use text question
    %fprintf('\n');
    for iLine = 1:length(question)-1
      fprintf([question{iLine} '\n']);
    end
    if toall
      % ask the question (with option for all)
      r = input([question{end} ' (y/n or a for Yes to all)? '],'s');
    else
      % ask question (without option for all)
      r = input([question{end} ' (y/n)? '],'s');
    end
  else
    if toall
      % verbose, use dialog
      r = questdlg(question,'','Yes','No','All','Yes');
      r = lower(r(1));
    else
      r = questdlg(question,'','Yes','No','Yes');
      r = lower(r(1));
    end
  end
  % make sure we got a valid answer
  if (lower(r) == 'n')
    r = 0;
  elseif (lower(r) == 'y')
    r = 1;
  elseif (lower(r) == 'a') & toall
    r = inf;
  else
    r =[];
  end
end


