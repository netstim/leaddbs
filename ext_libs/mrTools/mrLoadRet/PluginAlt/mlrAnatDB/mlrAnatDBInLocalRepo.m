% mlrAnatDBInLocalRepo.m
%
%        $Id:$ 
%      usage: tf = mlrAnatDBInLocalRepo(v)
%         by: justin gardner
%       date: 06/23/15
%    purpose: Tells whether the session is in the mlrAnatDB repo or not
%
%             v = newView;
%             tf = mlrAnatDBInLocalRepo(v);
%
function tf = mlrAnatDBInLocalRepo(v)

tf = false;

% get locations of session and repo
localRepoTop = mlrReplaceTilde(mrGetPref('mlrAnatDBLocalRepo'));
homeDir = mlrReplaceTilde(viewGet(v,'homeDir'));

% This is not a straight forward check since we 
% need to deal with links that can put directories
% in different places. So first find whether the
% directory of the repo top is in the homeDir
localRepoTopLoc = findstr(getLastDir(localRepoTop),homeDir);
if isempty(localRepoTopLoc),return,end

% it is, so now replace the path in homeDir with the localRepoTop
homeDirLocalRepo = fullfile(localRepoTop,homeDir(localRepoTopLoc+length(getLastDir(localRepoTop))+1:end));

% now check if the are the same directory by getting each ones link
% resolved path using pwd -P
curdir = pwd;
cd(homeDir);
[success homeDirResolved] = system('pwd -P');
cd(homeDirLocalRepo);
[success homeDirLocalRepoResolved] = system('pwd -P');
cd(curdir);

% see if they are the same or not
if ~strcmp(homeDirResolved,homeDirLocalRepoResolved)
  tf = false;
else
  tf = true;
end

