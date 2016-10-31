% mlrImageWriteNiftiHeader.m
%
%        $Id:$ 
%      usage: [hdr,fid] = mlrImageWriteNiftiHeader(hdr,fname,no_overwrite,leave_open)
%         by: justin gardner
%       date: 03/05/15
%    purpose: This is a transitional wrapper function that works identically
%             to cbiWriteNiftiHeader, so that we can transition away from 
%             cbiWriteNiftiHeader to mlrImage. You should
%             use mlrImageHeaderSave directly in most cases
%             and only use this function if you are 
%             replacing an old  call to cbiWriteNiftiHeader. This just
%             transforms inputs to way that mlrImage can handle and returns what
%             cbiWriteNiftiHeader would have. Ideally we would get rid of
%             cbiWriteNiftiHeader all together and work with mlrImage
%             structures which are meant to be file format independent
%             but I don't want to risk breaking everything... just yet.
%
%             Note that this does not support the leave_open
%             settings and does not return an fid and will warn if you
%             try to do any of that.
%
function [hdr,fid] = mlrImageWriteNiftiHeader(hdr,fname,no_overwrite,leave_open)

% check arguments
if ~any(nargin == [2:4])
  help mlrImageWriteNifti
  return
end

% return values
fid = [];

if nargout == 2
  disp(sprintf('(mlrImageWriteNiftiHeader) Returning open fid not supported. Returning empty fid'));
end

% check no_overwrite
if (nargin >= 3) && no_overwrite
  if isfile(fname) 
    disp(sprintf('(mlrImageWriteNiftiHeader) File %s exists. Not overwritng',fname));
    return
  end
end

% check leave_open
if (nargin >= 4) && no_overwrite
  disp(sprintf('(mlrImageWriteNiftiHeader) !!! leave_open not supported. Not returning an empty file handler'));
end

% call mlrImageHaderSave
mlrImageHeaderSave(fname,hdr);
