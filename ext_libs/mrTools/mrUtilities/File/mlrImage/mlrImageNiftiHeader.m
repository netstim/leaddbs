% mlrImageReadNifti.m
%
%        $Id:$ 
%      usage: function hdr = cbiReadNiftiHeader(fname)
%         by: justin gardner
%       date: 03/04/15
%    purpose: This is a transitional wrapper function that works identically
%             to cbiReadNiftiHeader, so that we can use mlrImage. 
%             Just transforms inputs to a way that mlrImage can
%             handle and return what cbiReadNiftiHeader would have. Ideally
%             we would get rid of cbiReadNiftiHeader all together and work with
%             mlrImage structures which are meant to be file format independent
%             but I don't want to break everything... just yet.
%
function hdr = mlrImageNiftiHeader(fname)

% check arguments
if ~any(nargin == [1])
  help mlrImageNiftiHeader
  return
end

% set this to check that mlrImage and cbiReadNifti return the
% same data
debugCheck = false;

% default to return empty
hdr = [];

% replace with a mlrImageLoad call
mlrHdr = mlrImageHeaderLoad(fname);
hdr = mlrImageGetNiftiHeader(mlrHdr);

% add some fields for filenames
if strcmp(mlrHdr.ext,'hdr')
  hdr.single_file = 0;
  hdr.hdr_name = setext(mlrHdr.filename,'hdr');
  hdr.img_name = setext(mlrHdr.filename,'img');
else
  hdr.single_file = 1;
  hdr.hdr_name = setext(mlrHdr.filename,'nii');
  hdr.img_name = setext(mlrHdr.filename,'nii');
end

% pixdim used to default to 1 for empty fields
hdr.pixdim(hdr.dim(1)+2:end) = 1;
hdr.dim(hdr.dim(1)+2:end) = 1;

% here is a check with old style cbiReadNiftiHeader call
if debugCheck
  % load old style
  hdr2 = cbiReadNiftiHeader(fname);
  % get fields
  fields = fieldnames(hdr2);
  theSame = true;
  for iField = 1:length(fields)
    if ~isfield(hdr,fields{iField})
      disp(sprintf('(mlrImageNiftiHeader) Field %s missing from mlrImage call',fields{iField}));
      theSame = false;
    elseif ~isequal(hdr.(fields{iField}),hdr2.(fields{iField}))
      disp(sprintf('(mlrImageNiftiHeader) Field %s different in mlrImage call',fields{iField}));
      theSame = false;
    end
  end
  % display match or not
  disp(sprintf('(loadTSeries) mlrImageNiftiHeader and cbiReadNiftiHeader return same data: %i',theSame));
  % incorrect then stop
  if ~theSame
    keyboard
  end
end


