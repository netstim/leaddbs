% mlrImageHeaderEdit.m
%
%        $Id:$ 
%      usage: h = mlrImageHeaderEdit(h)
%         by: justin gardner
%       date: 09/30/11
%    purpose: mrParamsDialog based gui for editing mlrImage headers
%
function h = mlrImageHeaderEdit(h)

% check arguments
if ~any(nargin == [1])
  help mlrImageHeaderEdit
  return
end

% make sure this is a header
[tf h] = mlrImageIsHeader(h);
if ~tf
  disp(sprintf('(mlrImageHeaderEdit) Must pass in an mlrImage header'));
  return
end

% if the matrices are empty fill them with nans so that 
% they show up correctly in the params dialog
if isempty(h.qform),qform = nan(4,4);else qform = h.qform;end
if isempty(h.sform),sform = nan(4,4);else sform = h.sform;end
if isempty(h.vol2mag),vol2mag = nan(4,4);else vol2mag = h.vol2mag;end
if isempty(h.vol2tal),vol2tal = nan(4,4);else vol2tal = h.vol2tal;end


% now gather details into paramsInfo
paramsInfo = {...
    {'filename',h.filename,'The filename'}...
    {'ext',h.ext,'The filename extension'}...
    {'type',h.type,'What type of file the image is (e.g. Nifti, data, fid, fdf, sdt, etc.)'}...
    {'nDim',h.nDim,'editable=0','Number of dimensions in the data'}...
    {'dim',h.dim,'editable=0','Size of data'}...
    {'pixdim',h.pixdim,'Size of data in mm. For volume dimension should be time in s'}...
    {'qform',qform,'Specifies the xform from the image coordinates to the magnet coordinates. Same as a nifti header qform when qform_code = 1. Will be empty if no such information exists'}...
    {'qformSet',{'Set qform','Identity','sform','vol2mag','vol2tal','clear'},'callback',@setXform,'callbackArg','qform','passParams=1','Set the qform'}...
    {'sform',sform,'Specifies the xform from the image coordinates to the magnet coordinates for the canonical volume. This is the same as a nifti header qform when sform_code = 1. Will be empty if no such info exists.'}...
    {'sformSet',{'Set sform','Identity','qform','vol2mag','vol2tal','clear'},'callback',@setXform,'callbackArg','sform','passParams=1','Set the sform'}...
    {'vol2mag',vol2mag,'This is the xform from the canonical volume coordinates to the magnet coordinates for the canonical. This is useful for checking which volume the image was aligned to (i.e. where it got the sform from). It is also useful for when trying to compute the talairach coordinates. Will be empty if no such info exists.'}...
    {'vol2magSet',{'Set vol2mag','Identity','qform','sform','vol2tal','clear'},'callback',@setXform,'callbackArg','vol2mag','passParams=1','Set the vol2mag'}...
    {'vol2tal',vol2tal,'This is the xform from the canonical volume coordinates to the talairach coordinates for the canonical. Note that we keep this matrix since it allows us to go back and change the talairach points on the canonical volume and then easily redo the tal xforma for the image by simply replacing this matrix. To compute the xform from the image coordinates to talairach coordiates, you can do: img2tal = vol2tal * inv(vol2mag) * sform * shiftOriginXform; Will be empty if no such info exists.'}...
    {'vol2talSet',{'Set vol2tal','Identity','qform','sform','vol2mag','clear'},'callback',@setXform,'callbackArg','vol2tal','passParams=1','Set the vol2tal'}...
	     };


% open up dialog box
params = mrParamsDialog(paramsInfo,'Edit mlrImage header');
if isempty(params),return,end

% keep the original, just in case the fields entered don't validate
originalHeader = h;

% change settings in header
h.filename = params.filename;
h.ext = params.ext;
h.type = params.type;
h.pixdim = params.pixdim;

% now do xform matrices 
if all(~isnan(params.qform(:)))
  h.qform = params.qform;
else
  h.qform = [];
end
if all(~isnan(params.sform(:)))
  h.sform = params.sform;
else
  h.sform = [];
end
if all(~isnan(params.vol2mag(:)))
  h.vol2mag = params.vol2mag;
else
  h.vol2mag = [];
end
if all(~isnan(params.vol2tal(:)))
  h.vol2tal = params.vol2tal;
else
  h.vol2tal = [];
end

% validate header
[tf h] = mlrImageIsHeader(h);
if ~tf
  disp(sprintf('(mlrImageHeaderEdit) Invalid header, aborting changes'));
  h = originalHeader;
  return
end

%%%%%%%%%%%%%%%%%%
%    setXform    %
%%%%%%%%%%%%%%%%%%
function retval = setXform(xformName,params)

retval = [];xform = [];
setTo = params.(sprintf('%sSet',xformName));

if strcmp(setTo,'Identity')
  % set to identity
  xform = eye(4);
  % see if we have pixdim
  if sum(~isnan(params.pixdim(1:3)))==3
    % then set diagonal to pixdim
    xform = diag([params.pixdim(1:3) 1]);
  end
elseif strcmp(setTo,'sform')
  xform = params.sform;
elseif strcmp(setTo,'vol2mag')
  xform = params.vol2mag;
elseif strcmp(setTo,'vol2tal')
  xform = params.vol2tal;
elseif strcmp(setTo,'qform')
  xform = params.qform;
elseif strcmp(setTo,'clear')
  xform = nan(4,4);
end


% set the xform
if ~isempty(xform)
  params.(xformName) = xform;
  % reset parameters
  params.(sprintf('%sSet',xformName)) = sprintf('Set %s',xformName);
  mrParamsSet(params);
end


