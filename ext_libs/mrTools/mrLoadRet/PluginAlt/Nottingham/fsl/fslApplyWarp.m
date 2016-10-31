% warpedData = fslApplyWarp(data, warpCoefFilename, tempInputFilename, tempRefFilename)
%
%   applies non-linear FSL registration to 3/4D data
%
% jb 23/07/2010
%
% $Id$
%
% input:
%      - data: 3/4D dataset, should be the same size as indicated in the header hdr
%      - warpCoefFilename
%      - tempFilename: name of the temporary NIFTI files what will get written on disk in order to call applywarp 
%      - hdr:  header giving the space of the data
%      - interpMethod
%      - verbose


function warpedData = fslApplyWarp(data, warpCoefFilename, tempFilename, hdr, interpMethod, verbose)


if ieNotDefined('interpMethod')
   interpMethod = mrGetPref('interpMethod');
end
if ieNotDefined('verbose')
   verbose = true;
end

cbiWriteNifti(tempFilename, data, hdr,[],[],[],verbose);


%convert matlab flags to FNIRT flags
switch(interpMethod)
  case 'nearest' %Nearest neighbor interpolation
    interpMethod = 'nn';
  case 'linear'%Linear interpolation (default)
    interpMethod = 'trilinear';
  case {'spline', 'cubic'}
    mrWarnDlg(['(fslApplyWarp) ' interpMethod ' not implemented in FNIRT, using ''sinc'''])
    interpMethod = 'sinc';
  case {'trilinear','nn','sinc'}
    %do nothing
  otherwise
    disp('(fslApplyWarp) Unknown interpolation method. Entering debug mode');
    keyboard;
    return;
end

try
  command =  sprintf('applywarp --ref=%s --in=%s --warp=%s --interp=%s --out=%s', tempFilename, tempFilename, warpCoefFilename, interpMethod, tempFilename);
  if verbose,disp(command);end;
  [s,w] = unix(command);

  if s ~=- 0 % unix error
    disp('UNIX error message:')
    disp(w)
    disp('-------------------')
    return
  end
catch 
  disp('(fslApplyWarp) There was a problem running FSL applywarp command')
  disp(sprintf('unix error code: %d; %s', s, w))
  return
end
%read the warped values
warpedData=mlrImageReadNifti(tempFilename);

if any(isnan(data(:))) %is there are NaNs, warp a mask of non-NaNs
  mask = ~isnan(data);
  warpedMask = logical(fslApplyWarp(mask, warpCoefFilename, tempFilename, hdr, 'nearest'));
  warpedData(~warpedMask) = NaN;
end

