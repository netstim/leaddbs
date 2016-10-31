% mlrImageWriteNifti.m
%
%        $Id:$ 
%      usage: [byteswritten,hdr]=mlrImageWriteNifti(fname,data,hdr,prec,subset,short_nan,verbose)
%         by: justin gardner
%       date: 03/05/15
%    purpose: This is a transitional wrapper function that works identically
%             to cbiWriteNifti, so that we can transition away from 
%             cbiWriteNifti to mlrImage. You should use mlrImageSave 
%             directly in most cases and only use this function if you are 
%             replacing an old  call to cbiWriteNifti. This Just transforms
%             inputs to way that mlrImage can handle and returns what
%             cbiWriteNifti would have. Ideally we would get rid of
%             cbiWriteNifti all together and work with mlrImage
%             structures which are meant to be file format independent
%             but I don't want to risk breaking everything... just yet.
%
%             Note that this does not support the prec, subset and short_nan
%             settings so will just print a warning and ignore them
%
function [byteswritten,hdr]=mlrImageWriteNifti(fname,data,hdr,prec,subset,short_nan,verbose)

% check arguments
if ~any(nargin == [1:7])
  help mlrImageWriteNifti
  return
end

% default arguments
if nargin < 2, data = [];end
if nargin < 3, hdr = [];end
if nargin < 4, prec = 'float32';end
if nargin < 5, subset = [];end
if nargin < 6, short_nan = [];end
if nargin < 7, verbose = false;end

% check that precision is set to float32 (since we don't deal
% with short data
if ~strcmp(lower(prec),'float32')
  disp(sprintf('(mlrImageWriteNifti) !!! Ignoring precision setting (%s not supported)',prec));
end

% check subset setting to be empty
if ~isempty(subset) && ~isequal(subset{:},[])
  disp(sprintf('(mlrImageWriteNifti) !!! Does not support saving an image subset. Saving all data'));
end

% check short_nan to be empty
if ~isempty(short_nan)
  disp(sprintf('(mlrImageWriteNifti) !!! Does not support short_nan setting. Ignoring.'));
end

% call mlrImageSave
mlrImageSave(fname,data,hdr);
  
  
  
  


