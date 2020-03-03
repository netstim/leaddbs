function [] = uw_overlay_bwthresh(anat_img, overlay_img, output_img, overlayThresh, overlayPercent, overlayAlpha, overlayInterp, overlayCoreg)

%   FUNCTION uw_overlay_bwthresh
%   Version 1.0
%
%   A function for merging overlay data with an anatomical image.  bwthresh
%   will do a binary threshold of the overlay data, and display it as a
%   single grayscale value on top of the anatomical.
%
%   DEPENDANCIES: spm5
%
%   SYNTAX:
%   [] = uw_overlay_bwthresh(anat_img, overlay_img, output_img,
%               overlayThresh, overlayPercent, overlayAlpha, overlayInterp, overlayCoreg)
%
%
%   *First three arguments are required, others optional.
%                  
%   anat_img - String: Anatomical image in ANALYZE or NIFTI format.
%              Assumes only a single data volume in image. Assumes that 
%              anat_img is the higher resolution image; the overlay_img will
%              be resampled to match its resolution.
%
%   overlay_img - String: Overlay image in ANALYZE or NIFTI format. Assumes 
%                 overlay is co-registered with anatomical, unless overlayCoreg
%                 is specified.
%
%   output_img - String: Specify filename for ANALYZE output image.
%
%   overlayThresh - Value to threshold the overlay_img.  Only results
%                   greater than the specified value will be overlaid
%                   Default is 0
%
%   overlayPercent - Specifies 'brightness' of overlay on the anatomical,
%                    as a percentage of the 'brightest' voxel found in the
%                    anatomical.  0 would be the value of the darkest
%                    voxel, 100 the brightest voxel.  You may go over 100%
%                    Default is 100
%
%   overlayAlpha - Specifies transparency of the overlay data.  100 is
%                  completely opaque, 0 is completely transparent.
%                  Default is 100
%
%   overlayInterp - Sets interpolation method for resampling the overlay
%                   data. Default is 1
%
%           0          Zero-order hold (nearest neighbour).
%           1          First-order hold (trilinear interpolation).
%           2->127     Higher order Lagrange (polynomial) interpolation using
%                      different holds (second-order upwards).
%          -127 - -1   Different orders of sinc interpolation.
%
%   overlayCoreg - A 4x4 affine rigid body transformation matrix that maps
%                  voxels in the overlay to voxels in the anatomy.  If not
%                  specified, coregistration is assumed,defaults to eye(4)
%                  Use spm_coreg or the coregistration interface in spm to
%                  determine this matrix.  Need to coregister the EPI to
%                  the anatomical, not the overlay data.
%
%   *Note: For most T1-weighted brain images the fat signal is much
%   brighter than gray or white matter, meaning larger than expected
%   values of overlayAlpha may be required to see the overlay.
%
%   Copyright 2007 Samuel A. Hurley - samuel.hurley[at]gmail.com
%   Universtiy of Wisconsin - Madison | Applied Neuro fMRI Lab
%   Published under GNU General Public License Version 3, see <http://www.gnu.org/licenses/>


%     This program is free software; you can redistribute it and/or modify
%     it under the terms of the GNU General Public License version 3 ONLY
%     as published by the Free Software Foundation
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.



% I. Error handling
% Test input arguments, check that the input files are valid
if nargin < 3  % Need 3 input args for filenames
    error('Arguments: You must specify at least three arguments.');
end

if nargin > 8 % Cannot have more than 8
    error('Arguments: You cannot specify more than eight arguments.  Check syntax.');
end


%Check that the first three arguments are strings & that the files exist:
if ~ ischar(anat_img)
    error('Arguments: anat_img must be a filename (string)');
end

if exist(anat_img) ~= 2     % exist returns 2 for a valid filename
    error(['Arguments:' anat_img ' - file not found.']);
end

if ~ ischar(overlay_img)
    error('Arguments: overlay_img must be a filename (string)');
end

if exist(overlay_img) ~= 2
    error(['Arguments:' overlay_img ' - file not found.']);
end

if ~ ischar(output_img)
    error('Arguments: output_img must be a filename (string)');
end

if exist(output_img) == 2  % If the output file already exists, we don't want to accidently overwrite
    error(['Arguments:' output_img ' - output file already exists.  Delete target file or change name of output.']);
end

% Check if optional arguments are supplied, otherwise use default values

switch nargin
    case 3      % No optional args, use defaults
        overlayThresh = 0;
        overlayPercent = 100;
        overlayAlpha = 100;
        overlayInterp = 1;
        overlayCoreg = eye(4);
    case 4
        overlayPercent = 100;
        overlayAlpha = 100;
        overlayInterp = 1;
        overlayCoreg = eye(4);
    case 5
        if overlayPercent < 0  % if the user sets a percentage below 0, default to 0
            overlayPercent = 0;
            disp('Warning: Overlay percentage below zero, defaulting to 0');
        end
        overlayAlpha = 100;
        overlayInterp = 1;
        overlayCoreg = eye(4);
    case 6
        if overlayPercent < 0  % if the user sets a percentage below 0, default to 0
            overlayPercent = 0;
            disp('Warning: Overlay percentage below zero, defaulting to 0');
        end
        if overlayAlpha > 100
            overlayAlpha = 100;
            disp('Warning: Overlay alpha set above 100.  Defaulting to 100');
        end
        if overlayAlpha < 1
            overlayAlpha = 1;
            disp('Warning: Overlay alpha set at or below zero.  Image will not be visible.');
        end
        overlayInterp = 1;
        overlayCoreg = eye(4);
    case 7
        if overlayPercent < 0  % if the user sets a percentage below 0, default to 0
            overlayPercent = 0;
            disp('Warning: Overlay percentage below zero, defaulting to 0');
        end
        if overlayAlpha > 100
            overlayAlpha = 100;
            disp('Warning: Overlay alpha set above 100.  Defaulting to 100');
        end
        if overlayAlpha < 1
            overlayAlpha = 0;
            disp('Warning: Overlay alpha set at or below zero.  Image will not be visible.');
        end
        overlayCoreg = eye(4);
    case 8
       % All arguments specified, just check for errors
        if overlayPercent < 0  % if the user sets a percentage below 0, default to 0
            overlayPercent = 0;
            disp('Warning: Overlay percentage below zero, defaulting to 0');
        end
        if overlayAlpha > 100
            overlayAlpha = 100;
            disp('Warning: Overlay alpha set above 100.  Defaulting to 100');
        end
        if overlayAlpha < 1
            overlayAlpha = 0;
            disp('Warning: Overlay alpha set at or below zero.  Image will not be visible.');
        end
end


% II: Declare variables, load stuff

% Load images with spm_vol
A = spm_vol(anat_img);      % A - anatomy image datastructure
B = spm_vol(overlay_img);   % B - overlay image datastructure

A = A(1,1);  % In case multiple volumes are stored, take only the first volume
B = B(1,1);  

C = spm_read_vols(A);       % C - anatomy image matrix

anatDim = A.dim;            % dimensions of anatomy image matrix
overlayDim = B.dim;         % dimensions of overlay image matrix
anatMat = A.mat;            % anatomy transformation matrix
overlayMat = B.mat;         % overlay transformation matrix

% Initialize other variables ahead of time
M = zeros(4);               % M - transformation matrix to scale up overlay onto anatomy
overlayScaled = zeros(anatDim);  % A scaled up version of overlay data
overlayMask   = zeros(anatDim);  % Binary mask of overlay data
anatMask      = zeros(anatDim);  % anat image, overlay elements excluded

D             = zeros(anatDim);  % D - the final overlay output image matrix


% III: Find maximum and minimum values in anat. image matrix C
elementsInC = anatDim(1)*anatDim(2)*anatDim(3);
Creshape    = reshape(C, 1, elementsInC);       % Use reshape to convert C into a vector

anatMaxval  = max(Creshape);  % anatMaxval - maximum value in anatomy image
anatMinval  = min(Creshape);  % anatMinVal - minimum value in anatomy image

% Determine the 'brightness' of the overlay, with specified overlayPercent
overlayPixelvalue = (overlayPercent/100) * (anatMaxval - anatMinval) + anatMinval;


% IV: Resample the overlay B to the same dimesions as the anat image C, using transform matrix M
voxelM = overlayMat \ overlayCoreg * anatMat;     % Transforms from voxel coords in anatomy to voxel coords in overlay
sliceM = 0;

for n = 1:anatDim(3)   % Loop through and reslice M
    sliceM(3) = n;
    M = voxelM * spm_matrix(sliceM);  % First shift to the nth slice in anat coords, then transform to overlay coords
    overlayScaled(:,:,n) = spm_slice_vol(B, M, anatDim(1:2), overlayInterp);
end

% V:  Calculate masks for anatomy and overlay image
overlayMask = overlayScaled >= overlayThresh; % Mask where overlay is greater than threshold
                                              % 0 for regions below thresh, 1 for regions at or
                                              % above threshold
anatMask = C .* (1 - overlayMask);  % Use overlay mask to remove elements from anat. image

% VI: Combine overlayMask and anatMask with overlayPixelvalue and overlayAlpha parameters
alphaBlend = round(C .* (1 - (overlayAlpha/100)) + overlayMask * overlayPixelvalue * (overlayAlpha/100));
D = anatMask + overlayMask .* alphaBlend;       % D is now our finished image matrix

% VII: Write to file
A.fname = output_img;
A.descrip = ['Overlay: threshold ' num2str(overlayThresh) ', ' num2str(overlayPercent) '% brightness, ' num2str(overlayAlpha) '% transparent.'];

% Don't use spm_write_vol - it will rescale the image data
% spm_write_vol(A, D);

A = spm_create_vol(A);  % Create new header
A = spm_write_plane(A, D, ':');  % Create new image
