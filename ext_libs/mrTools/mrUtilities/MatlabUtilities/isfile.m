% isfile.m
%
%      usage: isfile(filename)
%         by: justin gardner
%       date: 08/20/03
%       e.g.: isfile('filename')
%    purpose: function to check whether file exists
%
function [isit permission] = isfile(filename)

isit = 0;permission = [];
if (nargin ~= 1)
  help isfile;
  return
end
if isempty(filename),isit = 0;,return,end

% open file
fid = fopen(filename,'r');

% check to see if there was an error
if (fid ~= -1)
  % check to see if this is the correct path 
  if (~ispc && length(filename)>=1 && ~isequal(filename(1),filesep)) ... %checks whether filename was given as absolute or relative path
      || (ispc && length(filename)>=2 && (~isequal(filename(2),':') && ~ismember(filename(2),'/\'))) %on windows, an absolute path can include the disk identifier or use different separators 
      %if relative path, matlab might have opened a file with this name because it exists
      %somewhere in the path and not because it exists in the current directory (which is what we want to know)
%     openname = fopen(fid);
%     if ~strcmp(fullfile(pwd,filename),openname)
    if ~isempty(which(filename)) && ~strcmp(fullfile(pwd,filename),which(filename)) %using which instead of fopen, because older versions of matlab didn't return the full path
      %disp(sprintf('(isfile) Found file %s, but not in current path',openname));
      isit = false;
      fclose(fid);
      return;
    end
  end
  % close file and get permissions
  fclose(fid);
  [dummy permission] = fileattrib(filename);
  isit = 1;
else
  isit = 0;
end

