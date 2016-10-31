function volIC = intensityContrastCorrection(vol, crop, junkFrames);
% function volIC = intensityContrastCorrection(vol, crop);
%
% Intensity and contrast normalization. 
%
% Corrects for the intensity gradient applying a wiener-like filtering.
% Then removes mean intensity and normalizes by the contrast (standard
% deviation), both estimated by minimizing a robust error measure between
% the sampled histogram and a mixture of 2 gaussians. 
%
% If the input tseries is 4D, then intensity gradient, mean intensity and
% contrast are all estimated from the mean (over time) of the time series.
% Then the corrections are applied to each individual frame. This is the
% recommended usage rather than estimating the corrections individually for
% each frame.
%
% vol: input time series (x,y,z,t) or single volume (x,y,z).
%
% crop: optional crop region, 2x3 matrix as follows.
%          [startx starty startz;
%           endx   endy   endz]
%
% Based on mrAlign code originally written by Oscar Nestares and described
% in the Appendix of: Nestares O, Heeger DJ, Robust multiresolution
% alignment of MRI brain volumes. Magn Reson Med 43:705-715, 2000.
% Oscar Nestares 1999
% 1/2007 djh, updated to mrLoadRet-4.5
% $Id: intensityGradCorrection.m,v 1.1 2009/05/19 15:10:38 eli Exp $	


if ieNotDefined('crop')
    crop = [];
end
if ~isempty(crop)
    crop(1,:) = max(crop(1,:),[2 2 2]);
    crop(2,:) = min(crop(2,:),(size(vol(:,:,:,1)) - [2 2 2]));
end

if ieNotDefined('junkFrames'), junkFrames = 0; end
if ieNotDefined('doDriftCorr'),  doDriftCorr = 1; end
if ieNotDefined('doIntensityCorr'), doIntensityCorr = 1; end
if ieNotDefined('doContrastCorr'),  doContrastCorr = 1; end


% Check if 4D
nFrames = size(vol,4);

%%% drift correction
if doDriftCorr
  %%% Drift correction
  if (nFrames > 1) & doDriftCorr
    for f = 1:nFrames
      volMean(f) = mean(mean(mean( vol(crop(1,1):crop(2,1), crop(1,2):crop(2,2), crop(1,3):crop(2,3), f) )));
      vol(:,:,:,f) =   vol(:,:,:,f) / volMean(f);
    end
  end
end

% Initialize result to be a copy of input
volIC = vol;

if nFrames > 1
  vol = mean(vol(:,:,:,junkFrames+1:end), 4);
end

%%% Intensity correction
if doIntensityCorr
  % Estimate intensity gradient and noise
  [int noise] = estFilIntGrad(vol,0);

  % Crop noise and estimate noise variance
  if ~isempty(crop)
    noiseCrop = noise([crop(1,1):crop(2,1)],[crop(1,2):crop(2,2)],[crop(1,3):crop(2,3)]);
  else
    noiseCrop = noise;
  end
  II = find(~isnan(noiseCrop));
  sigma2 = mean(noiseCrop(II));
  
  % Correct intensity gradient applying a wiener-like filtering.
  % Old version:
  %    volC = vol.*int ./ (int.^2 + noise + sigma2/2);
  % Modified (djh, 1/2007) because noise is actually spatially homogeneous.
  denominator = (int.^2 + sigma2);
  volI = vol.*int ./ denominator;
  
  % Loop through frames, doing the same.
  for f=1:nFrames
    volIC(:,:,:,f) = (volIC(:,:,:,f) .* int) ./ denominator;
  end
end

%%% Contrast correction
if doContrastCorr
  % Crop intensity corrected volume
  if ~isempty(crop)
    volCrop = volI([crop(1,1):crop(2,1)],[crop(1,2):crop(2,2)],[crop(1,3):crop(2,3)]);
  else
    volCrop = volI;
  end
  
  % Build the histogram
  II = find(~isnan(volCrop));
  if isempty(II)
    disp(crop); 
    error('The cropped volume is empty');
  end
  % regHistogram is about 10x faster than hist
  [h x] = hist(volCrop(II), 256);
  % [h x] = regHistogram(volCrop(II), 256);
  % Normalize the histogram
  h = h/(sum(h)*mean(diff(x)));
  
  % Fit the histogram with a sum of two Gaussians using robust estimation.
  % Initial parameters: Choose the first gaussian with the actual mean and
  % variance, and the second gaussian around 1 and with 1/10 of the actual
  % variance, both with weights of 0.5 
  pinit = [mean(volCrop(II)) var(volCrop(II)) 1 var(volCrop(II))/10 0.5 0.5];
  % Minimize the robust measure of the error
  p = fminsearch('errGaussRob', pinit, [], x, h);
  
  % Select the mean closer to 1
  if abs(p(1)-1)>abs(p(3)-1)
    mu = p(3); sigma2 = p(4);
  else
    mu = p(1); sigma2 = p(2);
  end
  
  % Normalize, saturating for low and high values
  limit = 4;  % 4*std
  for f=1:nFrames
    frame = volIC(:,:,:,f);
    frame = (frame - mu)/sqrt(sigma2);
    frame = clip(frame,-limit,limit);
    volIC(:,:,:,f) = frame;
  end

end

return

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
if ~exist('lpf')
  lpf = conv([1 4 6 4 1]/16, [1 4 6 4 1]/16);
end

lpfZ = lpf;
if exist('PbyPflag')
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


function imb=addBorder(im, Nx, Ny, method);
% addBorder -  Adds a border to an image:
%	method=1 -> consider the image periodic
%	method=2 -> specular
%	method=3 -> repeat the edge pixel
%	
%	imb = addBorder(im,Nx,Ny,method);
%
% ON - 10/96 (from putborde)

[sy sx]=size(im);
imb=zeros(sy+2*Ny,sx+2*Nx);
imb(1+Ny:sy+Ny,1+Nx:sx+Nx)=im;

if method == 1
	imb(1:Ny,1+Nx:sx+Nx)=im(sy-Ny+1:sy,:);
	imb(sy+Ny+1:sy+2*Ny,1+Nx:sx+Nx)=im(1:Ny,:);
	imb(1+Ny:sy+Ny,1:Nx)=im(:,sx-Nx+1:sx);
	imb(1+Ny:sy+Ny,sx+Nx+1:sx+2*Nx)=im(:,1:Nx);
	imb(1:Ny,1:Nx)=im(sy-Ny+1:sy,sx-Nx+1:sx);
	imb(1:Ny,sx+Nx+1:sx+2*Nx)=im(sy-Ny+1:sy,1:Nx);
	imb(sy+Ny+1:sy+2*Ny,1:Nx)=im(1:Ny,sx-Nx+1:sx);
	imb(sy+Ny+1:sy+2*Ny,sx+Nx+1:sx+2*Nx)=im(1:Ny,1:Nx);
elseif method == 2
	imb(Ny:-1:1,1+Nx:sx+Nx)=im(2:Ny+1,:);
	imb(sy+2*Ny:-1:sy+Ny+1,1+Nx:sx+Nx)=im(sy-Ny:sy-1,:);
	imb(1+Ny:sy+Ny,Nx:-1:1)=im(:,2:Nx+1);
	imb(1+Ny:sy+Ny,sx+2*Nx:-1:sx+Nx+1)=im(:,sx-Nx:sx-1);
	imb(Ny:-1:1,Nx:-1:1)=im(2:Ny+1,2:Nx+1);
	imb(Ny:-1:1,sx+2*Nx:-1:sx+Nx+1)=im(2:Ny+1,sx-Nx:sx-1);
	imb(sy+2*Ny:-1:sy+Ny+1,Nx:-1:1)=im(sy-Ny:sy-1,2:Nx+1);
	imb(sy+2*Ny:-1:sy+Ny+1,sx+2*Nx:-1:sx+Nx+1)=im(sy-Ny:sy-1,sx-Nx:sx-1);
elseif method==3
	for k=1:Nx
		imb(Ny+1:sy+Ny,k)=im(:,1);
		imb(Ny+1:sy+Ny,k+sx+Nx)=im(:,sx);
	end
	for k=1:Ny
		imb(k,Nx+1:sx+Nx)=im(1,:);
		imb(k+sy+Ny,Nx+1:sx+Nx)=im(sy,:);
	end
else
	error('Not a valid value for method')
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Test/Debug for mrAlign
inp = ALIGN.inplanes;
inpIC = intensityContrastCorrection(inp);
inpICscaled = 32*(inpIC+4);
inpScaled = 256/(max(inp(:)) - min(inp(:))) * (inp - min(inp(:)));
slice = 10;
figure(1)
subplot(1,2,1)
image(squeeze(inpScaled(:,:,slice)));
colormap gray(256)
axis image
subplot(1,2,2)
image(squeeze(inpICscaled(:,:,slice)));
colormap gray(256)
axis image

% Test/Debug for mrLoadRet
view = newView;
scan = 1;
junkFrames = viewGet(view,'junkframes',scan);
nFrames = viewGet(view,'nframes',scan);
nSlices = viewGet(view,'nslices',scan);
tseries = loadTSeries(view,scan,'all');
tseries = tseries(:,:,:,junkFrames+1:junkFrames+nFrames);
meantseries = mean(tseries,4);

normMeantseries = intensityContrastCorrection(meantseries);
for s = 1:nSlices
    imagesc(normMeantseries(:,:,s));
    colormap(gray)
    axis image
    pause
end

normTseries = intensityContrastCorrection(tseries);
f = 1;
for s = 1:nSlices
    imagesc(normTseries(:,:,s,f));
    colormap(gray)
    axis image
    pause
end
s = round(nSlices/2);
for f = 1:nFrames
    imagesc(normTseries(:,:,s,f));
    colormap(gray)
    axis image
    drawnow
end

