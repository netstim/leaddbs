function im=ea_depth_of_field(fighandle)
if exist('fighandle','var')
    set(0,'CurrentFigure',fighandle);
    set(fighandle,'Visible','off');
end
camproj('persp');
im = [];
count = 0;
sam = 10;
w = 0.1;
target = camtarget;
mask = getframe(fighandle);

mask=double(mask.cdata(:,:,1));
mask(:)=1000;
mask(round(size(mask,1)/2),round(size(mask,2)/2))=-10000000000;
mask = imgaussfilt(mask,round(size(mask)/5));
mask=mask-min(mask(:));
mask=mask/max(mask(:));
offset=0.2;
mask=mask-offset;
mask(mask<0)=0;

mask=repmat(mask,1,1,3);

pos = campos;
[THETA,PHI] = meshgrid(w*linspace(-2*pi,2*pi,sam),w*linspace(-2*pi,2*pi,sam));
prev_theta = 0;
prev_phi = 0;
for s = randperm(numel(THETA))
    theta = THETA(s);
    phi = PHI(s);
    camorbit(-prev_theta,-prev_phi);
    prev_theta = theta;
    prev_phi = phi;
    camorbit(theta,phi);
    frame = getframe(fighandle);
    frame = im2double(frame.cdata);
    count = count+1;
    if isempty(im)
        im = frame;
    end
    im = im+frame.*mask;
end

%im = imtrim(im/count);
set(fighandle,'Visible','on');
im=im./((mask*count)+1);


function [out1,out2,out3,out4,out5] = imtrim(im,location)
%IMTRIM auto-crop an image like Photoshop's Edit>Trim feature, as of yet only
%grayscale imagse are supported
%
%   cropped = IMTRIM(IM) crop image based on top left corner
%
%   [cropped,t,b,l,r] = IMTRIM(IM) return cropped image and indices used to
%   crop the image. So cropped = im(t:b,l:r);
%
%   [t,b,l,r] = IMTRIM(IM) return only indices used to crop the image. So
%   cropped = im(t:b,l:r);
%
%   [...] = IMTRIM(IM,location) same as above but location may specify
%   top-left corner ('NorthWest') or bottom-right corner ('SouthEast') to be
%   the picel used in determining the auto-crop
%
%   Copyright Alec Jacobson, 2010
%

if(~exist('location'))
    location = 'NorthWest';
end

% gather corner value to which the image is compared
if(strcmp(location,'NorthWest'))
    corner_value = im(1,1);
elseif(strcmp(location,'SouthEast'))
    corner_value = im(1,1);
else
    error([location ' is not a valid location']);
end

% hard-coded threshold parameter, works equivalently with Photoshop's
% hardcoded parameter
threshold = 0.1;

% get difference of image with corner value
%difference = abs(im - corner_value)>0.1;
% should work for any number of channels
difference = sqrt(sum((im - corner_value).^2,3)) > ...
    sqrt(threshold^2*size(im,3));
[left_i,left_j] = ind2sub(size(difference),find(difference,1));
[right_i,right_j] = ind2sub(size(difference),find(difference,1,'last'));
[top_j,top_i] = ind2sub(size(difference'),find(difference',1));
[bottom_j,bottom_i] = ind2sub(size(difference'),find(difference',1,'last'));
if(nargout == 1)
    out1 = im(top_i:bottom_i,left_j:right_j);
elseif(nargout == 5)
    out1 = im(top_i:bottom_i,left_j:right_j);
    out2 = top_i;
    out3 = bottom_i;
    out4 = left_j;
    out5 = right_j;
else
    out1 = top_i;
    out2 = bottom_i;
    out3 = left_j;
    out4 = right_j;
end


function matrixOut = smooth2a(matrixIn,Nr,Nc)
% Smooths 2D array data.  Ignores NaN's.
%
%function matrixOut = smooth2a(matrixIn,Nr,Nc)
%
% This function smooths the data in matrixIn using a mean filter over a
% rectangle of size (2*Nr+1)-by-(2*Nc+1).  Basically, you end up replacing
% element "i" by the mean of the rectange centered on "i".  Any NaN
% elements are ignored in the averaging.  If element "i" is a NaN, then it
% will be preserved as NaN in the output.  At the edges of the matrix,
% where you cannot build a full rectangle, as much of the rectangle that
% fits on your matrix is used (similar to the default on Matlab's builtin
% function "smooth").
%
% "matrixIn": original matrix
% "Nr": number of points used to smooth rows
% "Nc": number of points to smooth columns.  If not specified, Nc = Nr.
%
% "matrixOut": smoothed version of original matrix
%
%
% 	Written by Greg Reeves, March 2009.
% 	Division of Biology
% 	Caltech
%
% 	Inspired by "smooth2", written by Kelly Hilands, October 2004
% 	Applied Research Laboratory
% 	Penn State University
%
% 	Developed from code written by Olof Liungman, 1997
% 	Dept. of Oceanography, Earth Sciences Centre
% 	G�teborg University, Sweden
% 	E-mail: olof.liungman@oce.gu.se

%
% Initial error statements and definitions
%
if nargin < 2, error('Not enough input arguments!'), end

N(1) = Nr;
if nargin < 3, N(2) = N(1); else N(2) = Nc; end

if length(N(1)) ~= 1, error('Nr must be a scalar!'), end
if length(N(2)) ~= 1, error('Nc must be a scalar!'), end

%
% Building matrices that will compute running sums.  The left-matrix, eL,
% smooths along the rows.  The right-matrix, eR, smooths along the
% columns.  You end up replacing element "i" by the mean of a (2*Nr+1)-by-
% (2*Nc+1) rectangle centered on element "i".
%
[row,col] = size(matrixIn);
eL = spdiags(ones(row,2*N(1)+1),(-N(1):N(1)),row,row);
eR = spdiags(ones(col,2*N(2)+1),(-N(2):N(2)),col,col);

%
% Setting all "NaN" elements of "matrixIn" to zero so that these will not
% affect the summation.  (If this isn't done, any sum that includes a NaN
% will also become NaN.)
%
A = isnan(matrixIn);
matrixIn(A) = 0;

%
% For each element, we have to count how many non-NaN elements went into
% the sums.  This is so we can divide by that number to get a mean.  We use
% the same matrices to do this (ie, "eL" and "eR").
%
nrmlize = eL*(~A)*eR;
nrmlize(A) = NaN;

%
% Actually taking the mean.
%
matrixOut = eL*matrixIn*eR;
matrixOut = matrixOut./nrmlize;












