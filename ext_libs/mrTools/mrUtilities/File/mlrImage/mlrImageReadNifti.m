% mlrImageReadNifti.m
%
%        $Id:$ 
%      usage: [data hdr] = mlrImageReadNifti(fname,subset,prec,short_nan,verbose)
%         by: justin gardner
%       date: 03/04/15
%    purpose: This is a transitional wrapper function that works identically
%             to cbiReadNifti, so that we can transition away from cbiReadNifti
%             to mlrImage. You should use mlrImageLoad directly in most
%             cases and only use this function if you are replacing an
%             old  call to cbiReadNifti. This Just transforms inputs to
%             a way that mlrImage can handle and returns what cbiReadNifti
%             would have. Ideally we would get rid of cbiReadNifit all
%             together and work with mlrImage structures which are meant
%             to be file format independent but I don't want to risk breaking
%             everything... just yet.
%
function [tseries hdr] = mlrImageReadNifti(fname,subset,prec,short_nan,verbose)

% check arguments
if ~any(nargin == [1:5])
  help mlrImageReadNifti
  return
end

% set this to check that mlrImage and cbiReadNifti return the
% same data
debugCheck = false;

if nargin == 1
  [tseries hdr] = cbiReadNifti(fname);
elseif nargin == 2
  [tseries hdr] = cbiReadNifti(fname,subset);
elseif nargin == 3
  [tseries hdr] = cbiReadNifti(fname,subset,prec);
elseif nargin == 4 
  [tseries hdr] = cbiReadNifti(fname,subset,prec,short_nan);
elseif nargin == 4
  [tseries hdr] = cbiReadNifti(fname,subset,prec,short_nan,verbose);
end
return
  
% default to return empty
data = [];hdr = [];

% check subset arguments (which make it so that you can return 
% some subset of the image). Default to return everything.
xMin = [];xMax = [];
yMin = [];yMax = [];
zMin = [];zMax = [];
volNum = [];
if nargin >= 2
  if length(subset) >= 1
    % get x arguments
    if ~isempty(subset{1}) 
      xMin = min(subset{1});
      xMax = max(subset{1});
    end
  end
  if length(subset) >= 2
    % get y arguments
    if ~isempty(subset{2}) 
      yMin = min(subset{2});
      yMax = max(subset{2});
    end
  end
  if length(subset) >= 3
    % get z arguments
    if ~isempty(subset{3}) 
      zMin = min(subset{3});
      zMax = max(subset{3});
    end
  end
  if length(subset) >= 4
    % get volNum arguments
    volNum = subset{4};
  end
end

% default prec
if (nargin < 3),prec = [];end

% short_nan not handled
if (nargin >= 4) && ~short_nan
  disp(sprintf('(mlrImageReadNifti) !!!! This wrapper does not handle setting the short_nan setting of cbiReadNifti to false (could implement, but why? -jg) THis may return unexpected results !!!!'));
else
  short_nan = [];
end

% default verbose
if (nargin < 5),verbose = [];end

% replace with a mlrImageLoad call
[tseries,hdr] = mlrImageLoad(fname,'xMin',xMin,'xMax',xMax,'yMin',yMin,'yMax',yMax,'zMin',zMin,'zMax',zMax,'volNum',volNum,'precision',prec,'nifti',true);

% here is a check with old style cbiReadNifti call.
%
% Actually this finds a bug in the cbiReadNifti code, which
% returns a non 1 based and incorrect time series (with some
% starting 0's if you do things like:
% t = loadTSeries(v,1,[1 2],[2 4],1,1);
% which should load slice 1 and 2, volume 2 and 4,
% but cbiReadNifti returns 0 followed by volumes 1:3
% note that other combinations work for cbiReadNifti like
% getting volumes starting at 1
% t = loadTSeries(v,1,[1 2],[1 4],1,1);
% or getting a single voxel
% t = loadTSeries(v,1,11,[2 4],15,18);
% so, probably never occurred in practice. jg
if debugCheck
  if ~isequal(exist('subset'),1),subset = [];end
  if ~isequal(exist('prec'),1),prec = [];end
  if ~isequal(exist('short_nan'),1),short_nan = [];end
  if ~isequal(exist('verbose'),1),verbose = [];end
  % load old style
  [tseries2,hdr2] = cbiReadNifti(fname,subset,prec,short_nan,verbose);
  % display match or not
  disp(sprintf('(loadTSeries) mlrImageLoad and cbiReadNifti return same data: %i',isequal(tseries,tseries2)));
  % incorrect then load the full series and check what is going on
  if ~isequal(tseries,tseries2)
    tseries3 = mlrImageLoad(fname,'nifti',true);
    keyboard
  end
end


