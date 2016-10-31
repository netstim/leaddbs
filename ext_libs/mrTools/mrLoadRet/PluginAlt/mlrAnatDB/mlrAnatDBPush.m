% mlrAnatDBPush.m
%
%        $Id:$ 
%      usage: mlrAnatDBPush(subjectID)
%         by: justin gardner
%       date: 07/05/15
%    purpose: Push the subjectID repo. 
%
function retval = mlrAnatDBPush(subjectID)

% check arguments
if ~any(nargin == [1])
  help mlrAnatDBPush
  return
end

% get pushType
pushType = lower(mrGetPref('mlrAnatDBPushType'));
if isempty(pushType),pushType = 'Normal';end
if isequal(pushType,'none')
  disp(sprintf('(mlrAnatDBPush) Push is set to none. Changes will not be uploaded to central repo'));
  return
end

% get the subjectID
subjectID = mlrAnatDBSubjectID(subjectID);
if isempty(subjectID),return,end

% find out if the repos exist
[localRepo localRepoLargeFiles] = mlrAnatDBGetRepo(subjectID,'noPull=1');

% keep current path
curpwd = pwd;

% push them if they exist
if ~isempty(localRepo)
  cd(localRepo)
  if isequal(pushType,'background')
    disppercent(-inf,sprintf('(mlrAnatDBPush) Pushing repo %s in the background. You should be able to work immediately, but if you shutdown matlab before the push has finished it may fail (in which case you should run mlrAnatDBPush again.',localRepoLargeFiles));
    mysystem(sprintf('hg push --new-branch &'));
  else
    disppercent(-inf,sprintf('(mlrAnatDBPush) Pushing repo %s',localRepo));
    mysystem(sprintf('hg push --new-branch'));
  end
  cd(curpwd);
  disppercent(inf);
end

% push them if they exist
if ~isempty(localRepoLargeFiles)
  cd(localRepoLargeFiles)
  if isequal(pushType,'background')
    disppercent(-inf,sprintf('(mlrAnatDBPush) Pushing repo %s in the background. You should be able to work immediately, but if you shutdown matlab before the push has finished it may fail (in which case you should run mlrAnatDBPush again.',localRepoLargeFiles));
    mysystem(sprintf('hg push --new-branch &'));
  else
    disppercent(-inf,sprintf('(mlrAnatDBPush) Pushing repo %s. This may take a few minutes',localRepoLargeFiles));
    mysystem(sprintf('hg push --new-branch'));
  end
  cd(curpwd);
  disppercent(inf);
end


%%%%%%%%%%%%%%%%%%
%    mysystem    %
%%%%%%%%%%%%%%%%%%
function [status,result] = mysystem(command)

disp(sprintf('(mlrAnatDBPut): %s',command));
[status,result] = system(command,'-echo');

