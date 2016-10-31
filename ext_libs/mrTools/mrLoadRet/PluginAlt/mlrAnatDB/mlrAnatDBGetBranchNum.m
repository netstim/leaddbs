% mlrAnatDBGetBranchNum.m
%
%        $Id:$ 
%      usage: branchNum = mlrAnatDBGetBranchNum(localRepo)
%         by: justin gardner
%       date: 06/24/15
%    purpose: Gets the branch number for the repository
%
%       e.g.: 
% localRepo = mlrAnatDBGetRepo(25);
% branchNum = mlrAnatDBBranchNum(localRepo);
%
function branchNum = mlrAnatDBGetBranchNum(localRepo)

branchNum = [];

% set path
curpwd = pwd;
cd(localRepo);

% get the branch name
[status,result] = mysystem(sprintf('hg branch'));
branchNameLoc = regexp(result,'v\d');
if ~isempty(branchNameLoc)
  branchNum = str2num(result(branchNameLoc+1:end));
else
  mrWarnDlg(sprintf('(mlrAnatDbGetBranchNum) Could not figure out version number. This should be the current branch of the repository and should be in the format vXXXX where XXX is a number (e.g. v0001). You can fix by going to repo %s and assiging a valid version number as the branch name'));
  cd(curpwd);
  return
end

cd(curpwd);

%%%%%%%%%%%%%%%%%%
%    mysystem    %
%%%%%%%%%%%%%%%%%%
function [status,result] = mysystem(command)

disp(sprintf('(mlrAnatDBGetBranchNum): %s',command));
[status,result] = system(command,'-echo');
