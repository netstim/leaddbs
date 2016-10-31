% mlrFixPermissions.m
%
%        $Id$ 
%      usage: mlrFixPermissions(dirname)
%         by: justin gardner
%       date: 10/15/09
%    purpose: Fix the permissions of an MLR directory
%
function mlrFixPermissions(dirname)

% check arguments
if ~any(nargin == [1])
  help mlrFixPermissions
  return
end

dirPermission = '775';
filePermission = '664';

if ~isfile(fullfile(dirname,'mrSession.mat'))
  if ~askuser(sprintf('(mlrFixPermissions) Directory %s is not an MLR session, continue anyway?',dirname))
    return
  end
end

% call the recursive function
disp('(mlrFixPermissions) Examining permissions...');
needsFix = recursivelyFixDirPermissions('',dirname,dirPermission,filePermission,1);

if needsFix & askuser('(mlrFixPermissions) Do you really want to change these permissions?')
  recursivelyFixDirPermissions('',dirname,dirPermission,filePermission,0);
elseif ~needsFix
  disp(sprintf('(mlrFixPermissions) Permissions for directory %s do not need fixing',dirname));
end
  

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   fixDirPermissions   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
function needsFix = recursivelyFixDirPermissions(dirpath,dirname,dirPermission,filePermission,dryrun)

needsFix = 0;

% get the current permissions for this directory
thisPermission = getPermissions(fullfile(dirpath,dirname));
% if permissions need to be changed
if ~isequal(thisPermission,dirPermission)
  if dryrun
    % just display what needs to be changed
    disp(sprintf('%s %s->%s',fullfile(dirpath,dirname),thisPermission,dirPermission));
  else
    % change the directory permission
    system(sprintf('chmod %s %s',dirPermission,fullfile(dirpath,dirname)));
  end
  needsFix = 1;
end

% now go through everything in the directory
dirlist = dir(fullfile(dirpath,dirname));
for i = 1:length(dirlist)
  if dirlist(i).isdir && ~isequal(dirlist(i).name,'.') && ~isequal(dirlist(i).name,'..')
    % recursively change everything in the directory
    needsFix = needsFix | recursivelyFixDirPermissions(fullfile(dirpath,dirname),dirlist(i).name,dirPermission,filePermission,dryrun);
  elseif ~dirlist(i).isdir
    % get the current permissions
    thisPermission = getPermissions(fullfile(dirpath,dirname,dirlist(i).name));
    % if permissions need to be changed
    if ~isequal(thisPermission,filePermission)
      if dryrun
	% just display what needs to be changed
	disp(sprintf('%s %s->%s',fullfile(dirpath,dirname,dirlist(i).name),thisPermission,filePermission));
      else
	% change the file permissions
	system(sprintf('chmod %s %s',filePermission,fullfile(dirpath,dirname,dirlist(i).name)));
      end
      needsFix = 1;
    end
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%
%%   getPermissions   %%
%%%%%%%%%%%%%%%%%%%%%%%%
function permissions = getPermissions(filename)

permissions = nan;
userPermissions = 0;groupPermissions = 0;worldPermissions = 0;
if isdir(filename)
  [ret val] = system(sprintf('ls -ld %s',filename));
else
  [ret val] = system(sprintf('ls -l %s',filename));
end
if ret == 0
  if length(val)>=10
    if val(2) == 'r';userPermissions = 4;end
    if val(3) == 'w';userPermissions = userPermissions+2;end
    if val(4) == 'x';userPermissions = userPermissions+1;end

    if val(5) == 'r';groupPermissions = 4;end
    if val(6) == 'w';groupPermissions = groupPermissions+2;end
    if val(7) == 'x';groupPermissions = groupPermissions+1;end

    if val(8) == 'r';worldPermissions = 4;end
    if val(9) == 'w';worldPermissions = worldPermissions+2;end
    if val(10) == 'x';worldPermissions = worldPermissions+1;end
    permissions = sprintf('%i%i%i',userPermissions,groupPermissions,worldPermissions);
  end
end
