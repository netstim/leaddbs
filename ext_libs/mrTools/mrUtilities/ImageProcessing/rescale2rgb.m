% rescale2rgb.m
%
%        $Id$
%      usage: rescale2rgb()
%         by: David Heeger
%       date: 10/16/07
%    purpose: taken out of refreshMLRDisplay
%
% function rgb = rescale2rgb(image,cmap,[clipMin,clipMax])
%
% Clips top and bottom tails of image histogram.
% Sets NaNs to the lowest index in the colormap
%
function rgb = rescale2rgb(image,cmap,clip,gamma)

% check arguments
if ~any(nargin == [2 3 4])
  help rescale2rgb
  return
end

if ~exist('clip','var')
  % Choose clipping based on histogram
  histThresh = length(image(~isnan(image)))/1000;
  [cnt, val] = hist(image(~isnan(image)),100);
  goodVals = find(cnt>histThresh);
  clipMin = val(min(goodVals));
  clipMax = val(max(goodVals));
  %if that doesn't work, take the 5th and 95th percentiles
  if clipMin==clipMax
    im = sort(image(~isnan(image)));
    clipMin = im(round(0.05*length(im)));
    clipMax = im(round(0.95*length(im)));
  end
  clip = [clipMin,clipMax];
else
  clipMin = clip(1);
  clipMax = clip(2);
end

% default is gamma of 1
if ~exist('gamma','var')
  gamma = 1;
end

% Clip
result = image;
result(image < clipMin) = clipMin;
result(image > clipMax) = clipMax;

% Scale
if isequal(clipMax-clipMin,0)
  rgb = clipMax*ones([size(image),3]);return
else
  indices = round((size(cmap,1)-1) * ((result-clipMin)/(clipMax-clipMin)).^gamma) + 1;
end
indices = max(1,min(indices,size(cmap,1)));

% Extract r,g,b components
r = zeros(size(image));
g = zeros(size(image));
b = zeros(size(image));
r(:) = cmap(indices,1);
g(:) = cmap(indices,2);
b(:) = cmap(indices,3);

% Stuff them into rgb
dims = [size(image),3];
rgb = cat(length(dims),r,g,b);


