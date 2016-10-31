function result = blur(im,levels,filt)
% BLUR: Blurs an image by blurring and subsampling repeatedly, followed by
% upsampling and blurring.
%
%      result=blur(im,[levels],[filt])
%
%      im - input image.
%      levels - number of times to blur and subsample (default is 1).
%      filt - blurring 1d filter to be applied separably to the rows and
%               cols of im (default ='binom5').
%
% DJH '96
% update 12/97 to conform to Eero's updated pyrTools

if ~exist('levels','var')
  levels=1;
end

if ~exist('filt','var')
  filt = 'binom5';
end

if isstr(filt)
  filt = namedFilter(filt);
end  

tmp = blurDn(im,levels,filt);

% save upBlurDEBUG tmp levels filt

result = upBlur(tmp,levels,filt);

% Make sure its the same size as the input image
result = result([1:size(im,1)],[1:size(im,2)]);

return;
