% mlrImageIsImage.m
%
%        $Id:$ 
%      usage: tf = mlrImageIsImage(filename)
%         by: justin gardner
%       date: 07/23/09
%    purpose: Test to see if the file is a valid image, if the extension is provided, will
%             test that filename explicitly. If no extension is given, will append the
%             mrGetPref('niftiFileExtension') to the filename
% 
%             returns 0 if not a valid image file, 1 otherwise
%
function tf = mlrImageIsImage(filename)

tf = false;

% check arguments
if ~any(nargin == [1])
  help mlrImageIsImage
  return
end

% put on nifti file extension, if there is no extension given
if isstr(filename) && isempty(getext(filename))
  filename = setext(filename,mrGetPref('niftiFileExtension'));
end

% no see if we can open the header
h = mlrImageHeaderLoad(filename);
if isempty(h),return,end

% check for empty image
if all(h.dim==0),return,end

% if we got here, it must be an image file
tf = true;
