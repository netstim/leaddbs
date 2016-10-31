% mlrImageHeaderSave.m
%
%      usage: mlrImageHeaderSave(filename,header)
%         by: justin gardner
%       date: 09/04/11
%    purpose: Saves the image header using filename. This
%             saves as a nifti header with a mat associated
%             file if there is info in the base
%
function [tf hdr] = mlrImageHeaderSave(filename,h)

% default return argument
tf = false;
hdr = [];

% check arguments
if ~any(nargin == [1 2])
  help mlrImageHeaderSave
  return
end

% set empty filename to filename in h
if (nargin < 2) || isempty(filename)
  filename = h.filename;
end

% check ext
ext = getext(filename);
if isempty(ext)
  filename = setext(filename,mrGetPref('niftiFileExtension'));
  ext = getext(filename);
end

% check for compressed file
compressFile = false;
if strcmp(ext,'gz')
  compressFile = true;
  % need to uncompress the file (so that we can just write the header)
  if isfile(filename)
    system(sprintf('gunzip %s',filename));
  end
  % strip off the gz
  filename = stripext(filename);
  % get the extension
  ext = getext(stripext(filename));
  % if the extension is empty then, we should set it to nii
  if isempty(ext)
    ext = 'nii';
    filename = setext(filename,ext);
  end
end

% make sure the header is valid
[tf h] = mlrImageIsHeader(h);
if ~tf
  disp(sprintf('(mlrImageHeaderSave) Invalid header'));
  return
end

% if we are not asked to pass back the header
% then save it
if nargout < 2,saveHeader = true;else,saveHeader = false;end

% no decide what type of header to make
switch(ext)
 case {'hdr','img','nii'}
  hdr = createHeaderNifti(h,filename,saveHeader);
 case {'sdt','spr','edt','epr'}
  hdr = createHeaderSDT(h,filename,saveHeader);
 otherwise
  disp(sprintf('(mlrImageHeaderSave) Unknown extension type: %s',ext));
  return
end

% now compress if asked for
if compressFile
  system(sprintf('gzip %s',filename));
end

% see if we need to save out a matlab extension
if ~isempty(h.base)
  % set fields that might have changed
  base = h.base;
  base.vol2mag = h.vol2mag;
  base.vol2tal = h.vol2tal;
  % and save
  matFilename = setext(filename,'mat');
  save(matFilename,'base');
end

% success
tf = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    createHeaderNifti    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
function hdr = createHeaderNifti(h,filename,saveHeader)

% create the nifti header
hdr = mlrImageGetNiftiHeader(h);

% save it
if saveHeader
  hdr = cbiWriteNiftiHeader(hdr,filename);
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%    createHeaderSDt    %
%%%%%%%%%%%%%%%%%%%%%%%%%
function hdr = createHeaderSDT(h,filename,saveHeader)

sdtFields = {'numDim',h.nDim,1;...
	     'dim',h.dim,1;...
	     'origin',zeros(size(h.dim)),0;...
	     'extent',h.pixdim(1:h.nDim).*(h.dim(h.nDim)/2-1)/10,0;...
	     'fov',h.pixdim(1:h.nDim).*h.dim(1:h.nDim)/10,0;...
	     'interval',h.pixdim,0;...
	     'dataType','REAL',0;...
	     'sdtOrient','NULL',0;...
	     'fidName','',0;...
	     'endian','ieee-le',1;...
	     'wordsize',4,0;...
	     'filelen',prod(h.dim)*4,0;...
	     'filename','',0;...
	    };

% copy fields form header if they exist
for i = 1:size(sdtFields,1)
  if isfield(h.hdr,sdtFields{i,1}) && ~sdtFields{i,3}
    % copy from header if it exists
    hdr.(sdtFields{i,1}) = h.hdr.(sdtFields{i,1});
  else
    hdr.(sdtFields{i,1}) = sdtFields{i,2};
  end
end

% add the orientation fields
if ~isempty(h.qform),hdr.qform = h.qform(:);end
if ~isempty(h.sform),hdr.sform = h.sform(:);end
if ~isempty(h.vol2mag),hdr.vol2mag = h.vol2mag(:);end
if ~isempty(h.vol2tal),hdr.vol2tal = h.vol2tal(:);end

% can't save the header by itself
if saveHeader
  disp(sprintf('(mlrImageHeaderSave) Cannot save header alone - must cal mlrImageSave'));
end

  