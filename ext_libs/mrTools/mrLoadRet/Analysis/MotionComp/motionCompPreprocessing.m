% motionCompPreprocessing.m
%
function [correctedTseries,crop,sliceTimes,baseVol,baseF] = motionCompPreprocessing(tseries,params,junkFrames,nFrames,totalFrames,sliceTimes,mask)

if ieNotDefined('mask'),mask = [];end

% % check arguments
% if ~any(nargin == [5])
%   help motionCompPreprocessing
%   return
% end

% initialize correctedTseries for preprocessing
correctedTseries = tseries;
baseFrame = params.baseFrame;
driftCorrection = params.driftCorrection;
gradIntensityCorrection = params.gradIntensityCorrection;
tSmooth = params.tSmooth;
crop = params.crop;
sliceTimeCorrection = params.sliceTimeCorrection;

% fix crop
if isempty(crop)
  crop = [2 2 2; size(tseries(:,:,:,1)) - [2 2 2]]; 
else
  crop(1,:) = max(crop(1,:),[2 2 2]);
  crop(2,:) = min(crop(2,:),(size(tseries(:,:,:,1)) - [2 2 2]));
end

% Get slice times and replicate the last frame of tseries for slice time
% correction
if sliceTimeCorrection
%  correctedTseries(:,:,:,end+1) = correctedTseries(:,:,:,end); %JB I commented this because the last frame should be replicated only after preprocessing (if it is at all replicated, see end of function) 
  switch params.sliceTimeString
    case 'end of TR'
      sliceTimes = sliceTimes - 1;  %JB: by removing 1TR to all the acquisition times, 
                                    %the last acquired frame becomes the one with the time closest to 0, hence the reference frame
                                    %(and the first-acquired slices will have to be interpolated closer to their respective next frame)
    case 'middle of TR'             
      sliceTimes = sliceTimes - 0.5; %JB: same with .5 TR
    case 'beginning of TR'         
      %JB: the first acquired frame is the reference time
      % so the slice acquisition times do not change (0 for the first acquired frame and close to 1 for the last-acquired frame,
      % which will have to be interpolated closer to the previous frame)
    otherwise
      mrErrorDlg('Invalid slice times');
  end
else
  sliceTimes = [];
end

% Drift correction
if driftCorrection
  for frame = 1:totalFrames
    frameMean(frame) = mean(mean(mean( correctedTseries(crop(1,1):crop(2,1), crop(1,2):crop(2,2), crop(1,3):crop(2,3), frame) )));
    correctedTseries(:,:,:,frame) =   correctedTseries(:,:,:,frame) / frameMean(frame);
  end
end

% correct for gradient in intensity over the volume
if gradIntensityCorrection
  tseriesMean = mean(correctedTseries(:,:,:,junkFrames+1:junkFrames+nFrames), 4);
  % Estimate intensity gradient and noise
  [int noise] = estFilIntGrad(tseriesMean,0);
  % Crop noise and estimate noise variance
  if ~isempty(crop)
    noiseCrop = noise([crop(1,1):crop(2,1)],[crop(1,2):crop(2,2)],[crop(1,3):crop(2,3)]);
  else
    noiseCrop = noise;
  end
  II = find(~isnan(noiseCrop));
  sigma2 = mean(noiseCrop(II));
  % Correct intensity gradient applying a wiener-like filtering.
  denominator = (int.^2 + sigma2);
  % Loop through frames, doing the same.
  for frame = 1:totalFrames
    correctedTseries(:,:,:,frame) = (correctedTseries(:,:,:,frame) .* int) ./ denominator;
  end
end

% temporal smoothing
if tSmooth
  % make a boxcar filter with tSmooth volumes on either side of each frame
  tmpTseries = correctedTseries;
  for frame = 1:totalFrames
    frameMin = frame - tSmooth;
    if frameMin < 1, frameMin = 1; end
    frameMax = frame + tSmooth;
    if frameMax > totalFrames, frameMax = totalFrames; end
    correctedTseries(:,:,:,frame) = nanmean(tmpTseries(:,:,:,frameMin:frameMax), 4);
  end
  clear tmpTseries;
end

% Get volume corresponding to base frame.  Other frames will be motion
% compensated to this one.
switch baseFrame
  case 'first'
    baseF = junkFrames+1;
    if sliceTimeCorrection && ~strcmp(params.sliceTimeString,'end of TR')
      baseF = max(2,baseF); % if slice-time correction, the  reference cannot be the first frame 
      %because that means the reference would be interpolated outside the time-series (unless everything is interpolated to the end of the TR)
    end
    baseVol = correctedTseries(:,:,:,baseF);
  case 'last'
    if sliceTimeCorrection && ~strcmp(params.sliceTimeString,'beginning of TR')
      baseF = size(correctedTseries,4)-1; % if slice-time correction, the  reference cannot be the last frame 
      %because that means the reference would be interpolated outside the time-series (unless everything is interpolated to the beginning of the TR)
    else
      baseF = size(correctedTseries,4);
    end
    baseVol = correctedTseries(:,:,:,baseF);
  case 'mean'
    baseF = 0;
    tseriesCrop = correctedTseries(:,:,:,junkFrames+1:junkFrames+nFrames);
    baseVol = nanmean(tseriesCrop,4);
  otherwise
    mrErrorDlg('Invalid base frame');
end

% apply mask if it exists
if ~isempty(mask) && any(mask(:)~=1)
  disp(sprintf('(motionCompPreprocessing) Applying mask'));
  % make sure mask size matches tSeries size
  if isequal(size(mask),size(correctedTseries(:,:,:,1)))
    % mask baseVol
    baseVol = mask.*baseVol;
    % and correctedTSeries
    for i = 1:size(correctedTseries,4)
      correctedTseries(:,:,:,i) = mask.*correctedTseries(:,:,:,i);
    end
  else
    disp(sprintf('(motionCompPreprocessing) Mask size %s does not match tSeries size %s',num2str(size(mask)),num2str(size(correctedTseries(:,:,:,1)))));
  end
end

%no need to replicate last slice as this is now done in warpaffineInterp3
% % %JB: replicate last slice only after all preprocessing
% % if sliceTimeCorrection
% %   correctedTseries(:,:,:,end+1) = correctedTseries(:,:,:,end);
% % end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Int, Noise] = estFilIntGrad(vol, PbyPflag, lpf);
% estFilIntGrad - Estimates the intensity gradient, using local mean
%
%    [Int, Noise] = estFilIntGrad(vol, <PbyPflag>, <lpf>);
%
% Inputs:
%  vol - input volume affected by the intensity gradient
%  PbyPflag - operates plane by plane if activated (default 0)
%  lpf - low pass filter (applied separably to x,y,z) used to
%        compute the local mean
%
% Outputs:
%  Int   - Estimated intensity
%  Noise - Spatial distribution of the noise, estimated as the
%          local variance
%
% Oscar Nestares - 5/99

% default low-pass filter
if ~exist('lpf', 'var')
  lpf = conv([1 4 6 4 1]/16, [1 4 6 4 1]/16);
end

lpfZ = lpf;
if exist('PbyPflag', 'var')
   if PbyPflag
      lpfZ = 1;
   end
end
B = (length(lpf)-1)/2;

% add border to the volume and estimate the intensity as the local mean
for k=1:size(vol,3)
	volB(:,:,k) = addBorder(vol(:,:,k), B, B, 2);
end
Int = convXYZsep(volB, lpf, lpf, lpfZ, 'repeat', 'valid');

% estimate the noise as the mean local variance
for k=1:size(Int,3)
   IntB(:,:,k) = addBorder(Int(:,:,k), B, B, 2);
end
Noise = convXYZsep((volB-IntB).^2, lpf, lpf, lpfZ, 'repeat', 'valid');

return
