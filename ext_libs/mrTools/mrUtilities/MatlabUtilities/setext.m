% setext.m
%
%        $Id$ 
%      usage: filename = setext(filename,ext,<override>)
%         by: justin gardner
%       date: 08/09/08
%    purpose: Makes sure that file the filename has the specified extension
%
%             If override is set to 1 (default) then any existing extension
%             will be changed to ext. If set to 0, then the extension
%             will be appended to the filename.
%
%             In either case, if the filename already has the extension
%             then filename will not be modified.
%
function filename = setext(filename,ext,override)

% check arguments
if ~any(nargin == [2 3])
  help setext
  return
end

% the default is to override whatever extension exists
if ieNotDefined('override'),override = 1;end

% strip any leading '.'
if length(ext) && isequal(ext(1),'.')
  ext = ext(2:end);
end

% set the extension
if override
  filename = sprintf('%s.%s',stripext(filename),ext);
else
  % add the extension to the filename, if the
  % extension isn't existing already
  if ~strcmp(getext(filename),ext)
    filename = sprintf('%s.%s',filename,ext);
  end
end

