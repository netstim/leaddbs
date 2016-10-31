% mlrImageSave.m
%
%      usage: mlrImageSave(filename,data,header)
%         by: justin gardner
%       date: 09/04/11
%    purpose: Saves the image. By default this will save a nifti
%             file. Will also save a associated mat file if there
%             is a base structure which contains info. These files
%             can be loaded with mlrImageLoad
%
function retval = mlrImageSave(filename,data,h)

% check arguments
if ~any(nargin == [2 3])
  help mlrImageSave
  return
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

%if data is a structure than grab its data field
if isstruct(data) 
  if isfield(data,'data')
    data = data.data;
  else
    disp(sprintf('(mlrImageSave) Unknown structure type passed as data'));
    return
  end
end

% create a header
if nargin < 3
  h = mlrImageHeaderLoad(data);
end

% make sure the header is valid
[tf h] = mlrImageIsHeader(h);
if ~tf
  disp(sprintf('(mlrImageSave) Invalid header'));
  return
end

% validate the dimensions
h.dim = size(data)';
h.nDim = length(h.dim);

% pass through mlrImageHeaderSave to get a nifti header
[tf hdr] = mlrImageHeaderSave(filename,h);
if ~tf
  disp(sprintf('(mlrImageSave) Could not make a valid header for saving'));
  return
end

switch (ext)
 case {'hdr','img','nii'}
  % write out nifti file
  cbiWriteNifti(filename,data,hdr);
 case {'sdt','spr','edt','epr'}
  hdr.data = data;
  writesdt(filename,hdr);
 otherwise
  disp(sprintf('(mlrImageSave) Unknown extension type: %s',ext));
end

% now compress if asked for
if compressFile
  system(sprintf('gzip %s',filename));
end

