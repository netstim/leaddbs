% getRelativePath.m
%
%        $Id$
%      usage: path2 = getRelativePath(path1,path2)
%         by: justin gardner
%       date: 12/18/07
%    purpose: gets path2 as a directory relative to path1
%             e.g.
%             getRelativePath('/x/y/z','/x/y/b') 
%             returns
%             '../b'
%
function newpath = getRelativePath(path1,path2)

% check arguments
if ~any(nargin == [2])
  help getRelativePath
  return
end

% remove any trailing separators from each
path1 = stripfilesep(path1);
path2 = stripfilesep(path2);

% go through and look for when the path 
% is different between path1 and path2
dir1 = ' ';dir2=' ';
while ~isempty(dir1) && strcmp(dir1,dir2)
  [dir1,path1] = strtok(path1,filesep);
  [dir2,path2] = strtok(path2,filesep);
end
if ~isempty(dir1)
  path1 = fullfile(filesep,dir1,path1);
end
if ~isempty(dir2)
  path2 = fullfile(filesep,dir2,path2);
end

% now count how many directories up in path1 you have to go
newpath = '';
[dir1,path1] = strtok(path1,filesep);
while ~isempty(dir1)
  [dir1,path1] = strtok(path1,filesep);
  newpath = fullfile(newpath,'..');
end

% then add what remains of the second path
newpath = fullfile(newpath,path2);

% remove file extension at beginning
newpath = fliplr(stripfilesep(fliplr(newpath)));