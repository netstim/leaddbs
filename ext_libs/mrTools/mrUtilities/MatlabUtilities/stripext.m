% stripext.m
%
%      usage: stripext(filename,<delimiter>)
%         by: justin gardner
%       date: 02/08/05
%    purpose: remove extension if it exists. if
%             delimiter is unspecified, uses .
%             i.e. returns 'filename' from 'filename.ext'
%
function retval = stripext(filename,delimiter)

if ~any(nargin == [1 2])
  help stripext;
  return
end

% dot delimits ext
if exist('delimiter')~=1,delimiter='.';,end

% if cell array, then run on each element of cell array
if iscell(filename)
  for iFilename = 1:length(filename)
    filename{iFilename} = stripext(filename{iFilename},delimiter);
  end
  retval = filename;
  return
end

% look for delimiter
retval = filename;
dotloc = findstr(filename,delimiter);

% if it is found then remove the last portion
% note that we could have used strtok, but that
% doesn't work for something stranget like file.name.ext
% it will return file instead of file.name 
if length(dotloc) > 0
  % make sure it does not occur before the last
  % file separator (which would occur if a directory
  % has the dot in it (i.e. hidden directory)
  fileseploc = findstr(filename,filesep);
  if isempty(fileseploc) || (dotloc(end)>fileseploc(end))
    retval = filename(1:dotloc(length(dotloc))-1);
  end
end
