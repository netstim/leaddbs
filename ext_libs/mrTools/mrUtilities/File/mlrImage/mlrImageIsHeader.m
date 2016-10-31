% mlrImageIsHeader.m
%
%        $Id:$ 
%      usage: [tf hdr] mlrImageIsHeader(h)
%         by: justin gardner
%       date: 08/19/11
%    purpose: Test whether the passed in strucutre is an
%              mlrImageHeader or not.
%
%              Basic fields for a mlrImage header are:
%
%              nDim: number of dimensions in image
%              dim: a row array with the dimensions of the image
%              pixdim: a row array with the dimensions in mm
%              qform: a 4x4 homogenous transform matrix that
%                     specifies the xform from the image
%                     coordinates to the magnet coordinates. Same
%                     as a nifti header qform when qform_code = 1.
%                     Will be empty if no such info exists.
%              sform: a 4x4 homogenous transform matrix that
%                     specifies the xform from the image
%                     coordinates to the magnet coordinates for
%                     the canonical volume. This is the same
%                     as a nifti header qform when sform_code = 1.
%                     Will be empty if no such info exists.
%            vol2mag: This is the xform from the *canonical* volume
%                     coordinates to the magnet coordinates for the
%                     canonical. This is useful for checking which
%                     volume the image was aligned to (i.e. where
%                     it got the sform from). It is also useful for
%                     when trying to compute the talairach coordinates.
%                     Will be empty if no such info exists.
%            vol2tal: This is the xform from the *canonical* volume
%                     coordinates to the talairach coordinates for the
%                     canonical. Note that we keep this matrix
%                     since it allows us to go back and change the
%                     talairach points on the canonical volume and
%                     then easily redo the tal xforma for the image
%                     by simply replacing this matrix. To compute
%                     the xform from the image coordinates to
%                     talairach coordiates, you can do:
%                     img2tal = vol2tal * inv(vol2mag) * sform * shiftOriginXform;
%                     Will be empty if no such info exists.
%               ext: The file extension for which the image was saved
%          filename: The original filename from which the image was loaded
%              type: The type of the image that was loaded 
%              base: An mlr base structure if one is associated
%                    with the file
%           talInfo: If non-empty contains the fields
%                    (AC,PC,SAC,IAC,PPC,AAC,LAC,RAC) which contain
%                    the coordinates of these talairach landmark points
%         
%
function [tf h] = mlrImageIsHeader(h)

tf = false;
% check arguments
if ~any(nargin == [1])
  help mlrImageIsHeader
  return
end

if (nargout == 2)
  % Add optional fields and return true if the image header with optional fields is valid
  requiredFields = {'dim'};
  optionalFields = {'qform',[];
		    'sform',[];
		    'vol2mag',[];
		    'vol2tal',[];
		    'base',[];
		    'hdr',[];
		    'ext','';
		    'type','unknown';
		    'filename','';
		    'talInfo',[];
		   };
  % The fields nDim and pixdim will be gereated below if necessary
else
  % Return 0 if the image header structure is missing any fields required or
  % optional (since w/out changing the header structure it is invalid).
  requiredFields = {'dim','nDim','pixdim','qform','sform','vol2mag','vol2tal','base','hdr','ext','type','filename','talInfo'};
  optionalFields = {};
end

% Initialize return value
tf = true;
if ~isstruct(h)
  tf = false;
  return
end

% check if it is nifti header
if mlrImageIsNiftiHeader(h)
  % If we are being asked to convert, then
  % grab fields from nifti header and convert into mlrHeader
  hdr.nDim = h.dim(1);
  hdr.dim = h.dim(2:h.dim(1)+1);
  hdr.qform = h.qform44;
  hdr.pixdim = h.pixdim(2:h.dim(1)+1);
  hdr.sform = h.sform44;
  hdr.hdr = h;
  % pump back through this function to add other fields
  [tf h] = mlrImageIsHeader(hdr);
  return
end

% Check required fields
for f = 1:length(requiredFields)
  fieldName = requiredFields{f};
  if ~isfield(h,fieldName)
    tf = false;
    return;
  end
end

% set nDim
if ~isfield(h,'nDim')
  h.nDim = length(h.dim);
end

% set pixdim to same length as dim with all ones
if ~isfield(h,'pixdim')
  h.pixdim = nan(size(h.dim));
end
h.pixdim(end+1:h.nDim) = nan;
h.pixdim = h.pixdim(1:h.nDim);

% Optional fields and defaults
for f = 1:size(optionalFields,1)
  fieldName = optionalFields{f,1};
  default = optionalFields{f,2};
  if ~isfield(h,fieldName)  
    h.(fieldName) = default;
  end
end

% make sure we have row arrays
h.dim = h.dim(:)';
h.pixdim = h.pixdim(:)';

% always have the same sort order of fields
h = orderfields(h);


%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   mlrImageIsHeader   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%
function tf = mlrImageIsNiftiHeader(h)

tf = false;
% a bit arbitrary, but just seeing if it has some expected fields
if all(isfield(h,{'intent_ps','qform44','sform44','dim','pixdim'}))
  tf = true;
end
