% freesurfer_load_xfm.m
%
%        $Id$ 
%      usage: xform = freesurfer_load_xfm(filename,<fieldname>)
%         by: justin gardner
%       date: 12/03/09
%    purpose: Loads a xfm matrix from a freesurfer file. For example:
%       e.g.: xform = freesurfer_load_xfm('mri/transforms/talairach.xfm')
%
function xform = freesurfer_load_xfm(filename,fieldname)

% check arguments
if ~any(nargin == [1 2])
  help freesurfer_load_xfm
  return
end

%mni2caret = loadXformFromFile('../mypreborder.sh');
%base2mni = loadXformFromFile('../FS003/mri/transforms/talairach.xfm');

% fieldname to look for
if nargin < 2
  fieldname = 'Linear_Transform';
end

xform = [];
if ~isfile(filename)
  disp(sprintf('(loadXformFromFile) Could not find xform file %s',filename));
  return
end
fieldname = lower(fieldname);
% open the file
fid = fopen(filename);
fline = lower(fgets(fid));
found = 0;
% go through file line by line
while ~isequal(fline,-1)
  % look for fieldname
  if ~found
    s = strtok(fline,' =');
    if strcmp(s,fieldname),found = 1;end
  else
    % grab rows of matrix until we get three lines
    if found < 4
      xform(found,:) = str2num(fline);
      found  = found+1;
    end
  end
  % load next line
  fline = lower(fgets(fid));
end
fclose(fid);

% check xform size
if ~isequal(size(xform),[3 4]) && ~isequal(size(xform),[4 4])
  disp(sprintf('(loadXformFromFile) Xform is not 3x4'));
  return
end

% add last row
if size(xform,1) == 3
  xform(4,:) = [0 0 0 1];
end

