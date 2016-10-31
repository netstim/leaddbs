% mlrAnatDBSetBranchNum.m
%
%        $Id:$ 
%      usage: tf = mlrAnatDBSetBranchNum(localREpo,branchNum)
%         by: justin gardner
%       date: 06/24/15
%    purpose: Sets the branch number of a repo (used for tracking changes as a group). This
%             is meant to be just a marker to keep the state of the repositories all in 
%             one convenient place (so you can go back to a state of the repo in which
%             say the rois and the session it was based on are all in the same state). In
%             that sense it is more like a tag - but tags are only for a changeset not for
%             the whole repo status so we use branches. The default branch should always
%             be the most recent branch which is the way mercurial seems to be working.
%
%             localRepo = mlrAnatDBGetRepo(25);
%             currentBranchNum = mlrAnatDBGetBranchNum(localRepo);
%             mlrAnatDBSetBranchNum(localRepo,currentBranchNum+1);
% 
%           
%
function tf = mlrAnatDBSetBranchNum(localRepo,branchNum)

% default return value
tf = false;

%check arguments
if nargin < 2
  help mlrAnatDBSetBranchNum;
  return
end

% set path
curpwd = pwd;
cd(localRepo);

% update branch number
branchName = sprintf('v%04i',branchNum);
[status,result] = system(sprintf('hg branch %s',branchName));
if status == 0
  tf = true;
end

cd(curpwd);

%%%%%%%%%%%%%%%%%%
%    mysystem    %
%%%%%%%%%%%%%%%%%%
function [status,result] = mysystem(command)

disp(sprintf('(mlrAnatDBSetBranchNum): %s',command));
[status,result] = system(command,'-echo');
