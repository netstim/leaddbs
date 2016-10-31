% mlrImageHeaderLoad.m
%
%        $Id: mlrImageHeaderLoad.m 1268 2008-08-19 15:43:36Z justin $ 
%      usage: mlrImageHeaderLoad(filename)
%         by: justin gardner
%       date: 08/19/08
%    purpose: loads a mlr image. This can handle various image types
%             including nifti. 
%
%             You can load using a filename, view/scanNum/groupNum,
%             select from a dialog box in the current or canonical
%             directory, or from a struct -> see mlrImageParseArgs
%             for details
%
%             For details on the header see:  mlrImageIsHeader
function retval = mlrImageHeaderLoad(varargin)

retval = [];

% check arguments
if nargin < 1
  help mlrLoadHeader
  return
end

% parse arguments
[imageArgs otherArgs] = mlrImageParseArgs(varargin);
if length(imageArgs) == 0
  disp(sprintf('(mlrImageHeaderLoad) No images specified to load'));
  return
end
verbose = [];
getArgs(otherArgs,{'verbose=0'});

% number of images headers to load. Note that for
% a single image, then we just return the header
% for multiple images, we will return a cell array
% of headers
nImages = length(imageArgs);allHeaders = {};

% cycle though each image argument
for iImage = 1:nImages
  header = [];
  % see if we have a structure with the filename and other
  % arguments fields
  loadArgs = {};altArgs = {};
  if isfield(imageArgs{iImage},'filename')
    % grab loadArgs which are passed to load routines
    if isfield(imageArgs{iImage},'loadArgs')
      loadArgs = imageArgs{iImage}.loadArgs;
    end
    % grab altArgs which are used in this routine
    if isfield(imageArgs{iImage},'altArgs')
      altArgs = imageArgs{iImage}.altArgs;
    end
    % get filename
    imageArgs{iImage} = imageArgs{iImage}.filename;
  end
  % if passed in a sturcutre with the field data
  % then grab the data and convert the rest of the
  % fields to a private header
  if isstruct(imageArgs{iImage})
    if isfield(imageArgs{iImage},'h')
      % see if it has an mlrImage header
      if mlrImageIsHeader(imageArgs{iImage}.h)
	header = imageArgs{iImage}.h;
      else
	% see if this is a nifti header
	header = mlrImageHeaderLoadPassedInNiftiHeader(imageArgs{iImage}.h);
      end
    end
    % if we didn't get it from above, then build a header
    % based on the data field
    if isempty(header) && isfield(imageArgs{iImage},'data')
      header.type = 'data';
      header.ext = '';
      header.hdr = rmfield(imageArgs{iImage},'data');
      header.dim = size(imageArgs{iImage}.data);
      [tf header] = mlrImageIsHeader(header);
    end
    % set in all headers if necessary
    if nImages > 1
      [tf allHeaders{end+1}] = mlrImageIsHeader(header);
    end
    continue;
  elseif isstr(imageArgs{iImage})
    filename = imageArgs{iImage};
    %check for file
    if ~isfile(filename) && ~isdir(filename)
      disp(sprintf('(mlrImageHeaderLoad) Could not find file %s',filename));
      if nImages > 1,allHeaders{end+1} = [];end
      continue
    end

    % set inital fields of header
    header.filename = filename;
    header.ext = getext(filename);
    
    % actually load from file, here is where we handle different
    % file types.
    switch lower(getext(filename))
     case {'img'}
      if ~isdir(filename)
	header = mlrImageHeaderLoadNifti(filename,header);
      else
	header = mlrImageHeaderLoadFDF(filename,header,verbose);
      end
     case {'hdr','nii'}
      header = mlrImageHeaderLoadNifti(filename,header);
     case {'gz'}
      header = mlrImageHeaderLoadCompressedNifti(filename,header);
     case {'sdt','spr','edt','epr'}
      header = mlrImageHeaderLoadSDT(filename,header);
     case {'fid'}
      header = mlrImageHeaderLoadFid(filename,header,verbose,loadArgs);
     otherwise
      disp(sprintf('(mlrImageHeaderLoad) Unknown image header type: %s',getext(filename)));
      return
    end

    % if we have an empty header it means we couldn't load the file
    if isempty(header)
      if nImages > 1,allHeaders{end+1} = [];end
      continue
    end

    % load the associated matlab header
    header.base = loadAssociatedMatlabHeader(filename);

    % set the vol2mag field - get it from the base 
    % struture if it has not already been set
    if ~isfield(header,'vol2mag')
      if isfield(header.base,'vol2mag')
	header.vol2mag = header.base.vol2mag;
      else
	header.vol2mag = [];
      end
    end
    % set the vol2tal field - get it from the base 
    % struture if it has not already been set
    if ~isfield(header,'vol2tal')
      if isfield(header.base,'vol2tal')
	header.vol2tal = header.base.vol2tal;
      else
	header.vol2tal = [];
      end
    end

    % set the talInfo field - get it from the base 
    % struture if it has not already been set
    if ~isfield(header,'talInfo')
      if isfield(header.base,'talInfo')
	header.talInfo = header.base.talInfo;
      else
	header.talInfo = [];
      end
    end

    % now fill in any missing fields from header
    [tf header] = mlrImageIsHeader(header);
    
  else
    disp(sprintf('(mlrImageHeaderLoad) Unknown argument'));
  end

  % keep all headers for returning
  if nImages > 1,[tf allHeaders{end+1}] = mlrImageIsHeader(header);end
end

% return a cell array if multiple images headers are to be loaded
% or just the one header otherwise
if nImages > 1
  retval = allHeaders;
else
  [tf retval] = mlrImageIsHeader(header);
end
  
   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    mlrImageHeaderLoadFDF    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function header = mlrImageHeaderLoadFDF(filename,header,verbose)

% read the nifti header
[d nifti] = fdf2nifti(filename,verbose,true);
if isempty(nifti),header = [];return;end
  
% set some info
header.type = 'fdf';
header.hdr = nifti;

% set rest of fields based on the nifti header
header = setHeaderBasedOnNifti(header,nifti);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    mlrImageHeaderLoadFid    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function header = mlrImageHeaderLoadFid(filename,header,verbose,loadArgs)

% read the nifti header
nifti = fid2niftihdr(filename,verbose,'loadArgs',loadArgs);
if isempty(nifti),header = [];return;end
  
% set some info
header.type = 'fid';
header.hdr = nifti;

% set rest of fields based on the nifti header
header = setHeaderBasedOnNifti(header,nifti);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    mlrImageHeaderLoadSDT    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function header = mlrImageHeaderLoadSDT(filename,header)

% load sdt
d = readsdt(filename,0);
if isempty(d),return,end

% set fileds
if any(strcmp(getext(filename),{'sdt','spr'}))
  header.type = 'sdt';
else
  header.type = 'edt';
end
header.hdr = d;

% set fields
header.nDim = length(d.dim);
header.dim = d.dim;
if isfield(d,'interval')
  header.pixdim  = d.interval(1:3)*10;
elseif isfield(d,'fov')
  header.pixdim = d.fov(1:3)./d.dim(1:3)*10
end

% set qform
if isfield(d,'qform')
  header.qform = reshape(d.qform,4,4);
end

% set sform
if isfield(d,'sform')
  header.sform = reshape(d.sform,4,4);
end

% set vol2mag and vol2tal
if isfield(d,'vol2mag')
  header.vol2mag = reshape(d.vol2mag,4,4);
end
if isfield(d,'vol2tal')
  header.vol2tal = reshape(d.vol2tal,4,4);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    mlrImageHeaderLoadCompressedNifti    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function header = mlrImageHeaderLoadCompressedNifti(filename,header)

header = [];

% make sure this is actually a nifti file that has been compressed
% by checking the filename
uncompressedFilename = stripext(filename);
if ~any(strcmp(getext(uncompressedFilename),{'nii'}))
  disp(sprintf('(mlrImageHeaderLoadCompressedNifti) File %s does not appear to be a compressed nifti file (which should have extensions like: filename.nii.gz)',filename));
  return
end

uncompressedExists = false;
if ~isfile(uncompressedFilename)
  % uncompress the file first
  system(sprintf('gunzip -c %s > %s',filename,uncompressedFilename));
else
  % uncompressed file already exists, so no need to gunzip
  uncompressedExists = true;
end

% read the nifti header
nifti = cbiReadNiftiHeader(uncompressedFilename);

% remove uncompressed (but only if it wasn't preexistent)
if ~uncompressedExists
  system(sprintf('rm -f %s',uncompressedFilename));
end

% check that it was loaded properly
if isempty(nifti)
  disp(sprintf('(mlrImageHeaderLoad) Could not load header for %s',uncompressedFilename));
  header=[];
  return;
end

% set some info
header.type = 'nifti';
header.hdr = nifti;

% set rest of fields based on the nifti header
header = setHeaderBasedOnNifti(header,nifti);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    mlrImageHeaderLoadNifti    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function header = mlrImageHeaderLoadNifti(filename,header)

% read the nifti header
nifti = cbiReadNiftiHeader(filename);
if isempty(nifti)
  disp(sprintf('(mlrImageHeaderLoad) Could not load header for %s',filename));
  header=[];
  return;
end

% set some info
header.type = 'nifti';
header.hdr = nifti;

% set rest of fields based on the nifti header
header = setHeaderBasedOnNifti(header,nifti);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   mlrImageHeaderLoadPassedInNiftiHeader   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function header = mlrImageHeaderLoadPassedInNiftiHeader(h)

header = [];
% first check some fields
checkFields = {'qform_code','qform44','sform_code','sform44','dim','pixdim'};
isNifti = true;
for iField = 1:length(checkFields)
  if isfield(h,checkFields{iField})
    header.(checkFields{iField}) = h.(checkFields{iField});
  else
    isNifti = false;
  end
end

% if we found all of the fileds then it is a nfti header
if isNifti
  header.type = 'nifti';
  header.hdr = h;
  header = setHeaderBasedOnNifti(header,h);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   setHeaderBasedOnNifti   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function header = setHeaderBasedOnNifti(header,hdr);

% set qform
if hdr.qform_code == 1
  header.qform = hdr.qform44;
elseif hdr.qform_code ~= 0
  header.qform = hdr.qform44;
  disp(sprintf('(mlrImageHeaderLoad) Unrecognized qform_code: %i',hdr.qform_code));
end

% set vol2mag and vol2tal based on sform
if hdr.sform_code == 1
  header.sform = hdr.sform44;
elseif hdr.sform_code ~= 0
  header.sform = hdr.sform44;
  disp(sprintf('(mlrImageHeaderLoad) Unrecognized sform_code: %i',hdr.sform_code));
end  

% set dimensions
header.nDim = hdr.dim(1);

% check nDIm - should not be set to > 7
if header.nDim > 7
  header.nDim = find(hdr.dim==0);
  if isempty(header.nDim)
    header.nDim == 7;
  else
    header.nDim = header.nDim(1);
  end
  disp(sprintf('(mlrImageheaderLoad) Nifti header has the first element of dim set to %i - but the max number of dimensions possible is 7. Resetting to %i',hdr.dim(1),header.nDim));
end

% set the rest of the dimensions
header.dim = hdr.dim(2:header.nDim+1);
header.pixdim = hdr.pixdim(2:header.nDim+2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    loadAssociatedMatlabHeader    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function base = loadAssociatedMatlabHeader(filename)

base = [];
matlabFilename = setext(filename,'mat');
if isfile(matlabFilename)
  % load the header
  matHeader = load(matlabFilename);
  % check for field
  if ~isfield(matHeader,'base')
    disp(sprintf('(mlrImageHeaderLoad) Ignoring associated mat file %s because is not a MLR header',filename));
  else
    matHeader.base.data = [];
    [tf header] = isbase(matHeader.base);
    if ~tf
      disp(sprintf('(mlrImageHeaderLoad) Ignoring associated mat file %s because is not a MLR header',filename));
    else
      base = matHeader.base;
    end
  end
end
