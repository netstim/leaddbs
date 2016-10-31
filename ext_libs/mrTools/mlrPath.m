% mlrPath.m
%
%        $Id:$ 
%      usage: mlrPath
%         by: justin gardner
%       date: 09/09/14 (from nott.m)
%    purpose: allows happy disambiguation of paths for people
%             with mrTools and vista installed. Examines what
%             paths you have installed and will shut down
%             conflicting paths from vista so that you can run
%             mrTools - also, will bring back those paths when
%             you quit and go back
%
%             run without any arguments to switch between mrTools/vistasoft
%
%             mlrPath
%
%             Or set which one you want to use:
%
%             mlrPath mrTools
%             % -or-
%             mlrPath vistasoft
%
%             figure out which mrTools path you are using
% 
%             mlrPath which
%
%             Or switch back to the one you started with
%
%             mlrPath revert
%
%
function mlrPath(switchTo)

verbose = true;

% try to get paths for vista and for mlr
mlrRoot = getRootPath('mlr');
vistaRoot = getRootPath('vista');

% try to guess vista from mlr if we didn't get it
if isempty(vistaRoot) && ~isempty(mlrRoot)
  vistaRoot = fullfile(fileparts(mlrRoot),'vistasoft');
  if ~isdir(vistaRoot), vistaRoot = []; end
end

% try to guess mlr from vista if we didn't get it
if isempty(mlrRoot) && ~isempty(vistaRoot)
  mlrRoot = fullfile(fileparts(vistaRoot),'mrTools');
  if ~isdir(mlrRoot), mlrRoot = []; end
end

% add matlab utilities if we have it
matlabUtilitiesPath = fullfile(mlrRoot,'mrUtilities');
if ~isempty(mlrRoot) && isempty(findstr(matlabUtilitiesPath,path))
  addpath(genpath(matlabUtilitiesPath));
end

% save as prefs
if ~exist('mrSetPref')==2
  if ~isempty(mlrRoot), mrSetPref('mlrPath',mlrRoot);end
  if ~isempty(vistaRoot), mrSetPref('vistaPath',vistaRoot);end
end

% get which one we are using
whichMLR = fileparts(fileparts(which('mrAlign')));
[dump whichMLRName] = fileparts(whichMLR);
whichMLRName = lower(whichMLRName);

% check for secondary path
hasBothPaths = false;
otherPath = '';
if strcmp(whichMLR,mlrRoot)
  if ~isempty(which('mrVista')),hasBothPaths = true;end
  otherPath = vistaRoot;
else
  if ~isempty(which('mrLoadRet')),hasBothPaths = true;end
  otherPath = mlrRoot;
end
  
% decide which one to switch to.
if (nargin <= 0) || isempty(switchTo)
  if strcmp(whichMLR,vistaRoot) 
    switchTo = mlrRoot;
  elseif strcmp(whichMLR,mlrRoot) 
    switchTo = vistaRoot;
  end
end

% get the top level names of mlrRoot
if ~isempty(mlrRoot)
  [dump mlrName] = fileparts(mlrRoot);
  if isempty(mlrName),mlrName = mlrRoot;end
  mlrName = lower(mlrName);
else
  % set root to empty string (so it works in fileparts)
  mlrName = '';mlrRoot = '';
end
% get the top level names of vistaRoot
if ~isempty(vistaRoot)
  [dump vistaName] = fileparts(vistaRoot);
  if isempty(vistaName),vistaName = vistaRoot;end
  vistaName = lower(vistaName);
else
  vistaName = '';vistaRoot = '';
end
if ~isempty(switchTo)
  [dump switchToName] = fileparts(switchTo);
  if isempty(switchToName),switchToName = switchTo;end
  switchToName = lower(switchToName);
else
  switchToName = '';switchTo = '';
end

% a few names that the user can put for switchTo
% make these the same as the actual directory name
if any(strcmp(switchToName,{'mlr','mrtools'}))
  switchToName = mlrName;
  switchTo = mlrRoot;
end
if any(strcmp(switchToName,{'vista','vistasoft'}))
  switchToName = vistaName;
  switchTo = vistaRoot;
end

% user wants to revert to original path that was set
% this happens when mrTools quits
addOtherPath = [];
if any(strcmp(switchToName,{'revert'}))
  global mlrOriginalPath;
  if isempty(mlrOriginalPath)
    % not set in global, so grab from lastPath preference
    revertPath = mrGetPref('lastPath');
  else
    revertPath = mlrOriginalPath;
  end
  if isempty(revertPath)
    return
  end
  % get the path parameters
  [switchTo switchToName] = fileparts(revertPath);
  switchTo = fullfile(switchTo,switchToName);
  switchToName = lower(switchToName);
  % check whether we have to install the other path as well
  global mlrOriginalHasBothPaths;
  if isequal(mlrOriginalHasBothPaths,true)
    if isequal(switchToName,mlrName)
      addOtherPath = vistaRoot;
    else
      addOtherPath = mlrRoot;
    end
  end
end
  
% unknown switchTo or missing path
if ~any(strcmp(switchToName,{mlrName,vistaName}))
  disp(sprintf('(mlrPath) Currently using: %s',whichMLR));
  % check both paths
  if hasBothPaths
    disp(sprintf('(mlrPath) Secondary path is: %s',otherPath));
  end
  return
end  

% no vista path
if (strcmp(switchToName,vistaName) && isempty(vistaRoot))
  disp(sprintf('(mlrPath) Could not find vistasoft path'));
  return
end

% already have correct path (and there is only one path)
% so nothing to do
if isequal(whichMLRName,switchToName) && ~hasBothPaths && isempty(addOtherPath)
  return
end

% display what we are going to do
if verbose
  % check both paths
  if hasBothPaths
    disp(sprintf('(mlrPath) Had both %s and %s leaving only %s',whichMLR,otherPath,switchTo));
  else
    disp(sprintf('(mlrPath) Switching from %s to %s',whichMLR,switchTo));
  end
end

% now save what path we are using
global mlrOriginalPath
global mlrOriginalHasBothPaths
if isempty(mlrOriginalPath)
  mlrOriginalPath = whichMLR;
  mlrOriginalHasBothPaths = hasBothPaths;
end
mrSetPref('lastPath',whichMLR);

% first remove everyone
warning('off','MATLAB:rmpath:DirNotFound')
rmpath(genpath(mlrRoot));
rmpath(genpath(vistaRoot));
warning('on','MATLAB:rmpath:DirNotFound')

% if adding the other path (this happens when we are adding
% back both paths since this is the way they used to be - i.e. 
% a revert call. Note this comes first so that it will not have
% precedence in the path
if ~isempty(addOtherPath)
  addpath(genpath(addOtherPath));
  if verbose
    disp(sprintf('(mlrPath) Adding secondary path: %s',addOtherPath));
  end
end

% switch path
if strcmp(switchToName,mlrName)
  % add mrTools
  addpath(genpath(mlrRoot));
  % selectively add some paths from vista
  pathNames = {'mrDiffusion','mrMesh','utilities','mrBOLD/Utilities','fileFilters','external/pyrTools'};
  addpath(vistaRoot);
  for i = 1:length(pathNames)
    thisPath = fullfile(vistaRoot,pathNames{i});
    if isdir(thisPath)
      addpath(genpath(thisPath));
    end
  end
elseif strcmp(switchToName,vistaName)
  % add vista
  addpath(genpath(vistaRoot));
  % add just the top level directory for mlr (so that we can have this
  % function)
  addpath(mlrRoot);
else
  % add back both
  addpath(genpath(vistaRoot));
  addpath(genpath(mlrRoot));
end

%%%%%%%%%%%%%%%%%%%%%
%    getRootPath    %
%%%%%%%%%%%%%%%%%%%%%
function rootPath = getRootPath(packageName)

% get the function that returns the path of the package and pref name
functionName = sprintf('%sRootPath',packageName);
prefName = sprintf('%sPath',packageName);

% see if that function is on the path
if exist(functionName)
  rootPath = eval(functionName);
else
  rootPath = [];
  % try to get the path
  if exist('mrGetPref')==2
    rootPath = mrGetPref(prefName);
  end
end

% make sure it exists
if ~isdir(rootPath),rootPath = [];end
