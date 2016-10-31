% mlrAnatDBPreferences.m
%
%        $Id:$ 
%      usage: mlrAnatDBPreferences()
%         by: justin gardner
%       date: 06/30/15
%    purpose: Bring up dialog box to set preference for mlrAnatDB. Note that this
%             is setup to take (but not use) arguments since it can be called from a callback
%
function mlrAnatDBPreferences(varargin)

% get repo locations
centralRepo = mrGetPref('mlrAnatDBCentralRepo');
localRepoTop = mrGetPref('mlrAnatDBLocalRepo');

% set defaults
if isempty(centralRepo),centralRepo = '';end
if isempty(localRepoTop),localRepoTop = '~/data/mlrAnatDB';end

% get lockLocal
lockLocal = mrGetPref('mlrAnatDBLockLocal');
if isempty(lockLocal)
  lockLocal = false;
  mrSetPref('mlrAnatDBLockLocal',lockLocal,false);
end

% get push type
pushType = mrGetPref('mlrAnatDBPushType');
pushTypes = {'Normal','Background','None'};
if isempty(pushType)
  pushType = 'Normal';
  mrSetPref('mlrAnatDBPushType',pushType,false);
end

% get wiki location
wikiHostname = mrGetPref('mlrAnatDBWikiHostname');
wikiDirname = mrGetPref('mlrAnatDBWikiDirname');

% setup params info for mrParamsDialog
paramsInfo = {...
    {'mlrAnatDBCentralRepo',centralRepo,'Location of central repo, Typically on a shared server with an https address (or via ssh), but could be on a shared drive in the file structure.'}...
    {'mlrAnatDBLocalRepo',localRepoTop,'Location of local repo which is typically under a data directory - this will have local copies of ROIs and other data but can be removed the file system as copies will be stored in the central repo'}...
    {'mlrAnatDBLockLocal',lockLocal,'type=checkbox','If this is clicked on then your local repository will be locked, meaning that it will not pull from the central database and thus no longer be updated. This is useful if you are working on a project at a point in which you do not want to change any ROIs'}...
    {'mlrAnatDBPushType',putOnTopOfList(pushType,pushTypes),'Sets how you want to push to the central repo. Normal will block execution until the push has finished (default and recommended behavior). Background will push as a background process (note that if you shutdown matlab or the shell before the push is completed, it will stop. None means to never push'}...
    {'mlrAnatDBWikiHostname',wikiHostname,'Hostname of dokuwiki wiki for logging commits. Leave blank if you do not have a wiki. '}...
    {'mlrAnatDBWikiDirname',wikiDirname,'Directory on mlrAnatDBWikiHostname where dokuwiki is for logging commits. Leave blank if you do not have a wiki. '}...
};

% and display the dialog
params = mrParamsDialog(paramsInfo);

% save params, if user did not hit cancel
if ~isempty(params)
  mrSetPref('mlrAnatDBCentralRepo',params.mlrAnatDBCentralRepo,false);
  mrSetPref('mlrAnatDBLocalRepo',params.mlrAnatDBLocalRepo,false);
  mrSetPref('mlrAnatDBLockLocal',params.mlrAnatDBLockLocal,false);
  mrSetPref('mlrAnatDBPushType',params.mlrAnatDBPushType,false);
  mrSetPref('mlrAnatDBWikiHostname',params.mlrAnatDBWikiHostname,false);
  mrSetPref('mlrAnatDBWikiDirname',params.mlrAnatDBWikiDirname,false);
end

