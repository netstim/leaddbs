function [balanced,cmap] = ea_autocontrast(img, percent)
% Balance the contrast for visualization of nifti image
%
% Follow the fashion in MRIcron, the brightest PERCENT% of the voxels are
% shown as the maximum intensity, and the darkest PERCENT% of the voxels
% are shown as the minimum intensity. A gray color map is also calculated 
% to be used in myslicer.

if nargin < 2
    percent = 5;
end

[lBound, uBound] = percentbound(img, percent);

balanced = single(img);

balanced(balanced>=uBound) = uBound;
balanced(balanced<=lBound) = lBound;

cmap = gray(round(uBound-lBound));


function [lBound, uBound] = percentbound(X, percent)
% Calculate the lower and upper bound of the input X
% uBound: the minimum value of the largest PERCENT% of the elements in X
% lBound: the maximum value of the smallest PERCENT% of the elements in X

sorted = sort(single(X(:)));

lBound = sorted(round(numel(sorted)*percent/100)+1);

uBound = sorted(round(numel(sorted)*(1-percent/100))-1);
