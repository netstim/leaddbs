function out = upSampleRep(im,resSize)
% out = upSampleRep(im, nLevels,resSize)
%
% Upsamples the image im by pixel replication.
%
% djh, 2/19/2001

upSampleFactor = resSize./size(im);
y = [0:resSize(1)-1];
x = [0:resSize(2)-1];
ysub = floor(y/upSampleFactor(1));
xsub = floor(x/upSampleFactor(2));
out = im(ysub+1,xsub+1);

return

%%% Debug/test

im = zeros(3,3);
im(2,2)=1;

upIm = upSampleRep(im,[6,6])
upIm = upSampleRep(im,[9,9])
upIm = upSampleRep(im,[6,7])

