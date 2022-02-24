function [overlap1, normOverlap1, overlap2, normOverlap2] = ea_nii_overlap(image1, image2, binary, threshold)
% Calculate the overlap between two images (in the same space)
%
% Return the overlaps in VOXEL (NOT IN MM^3), and normalized overlaps (
% normalized by the voxels of nii1 and nii2 respectively):
%     overlap1 and normOverlap1 are calculated using image1 as reference.
%     overlap2 and normOverlap2 are calculated using image2 as reference.

% Binarizing images by default
if ~exist('binary', 'var')
    binary = 1;
end

% Threshold set to nan by default
if ~exist('threshold', 'var')
    threshold = nan;
end

% Reslice image2 using image1 as reference
reslicedImage2 = strrep(image2, '.nii', '_resliced.nii');
ea_fsl_reslice(image2, image1, reslicedImage2, 'trilinear', 0);

nii1 = ea_load_nii(image1);
nii2 = ea_load_nii(reslicedImage2);

% Binarize when needed
if binary
    [nii1.img, nii2.img] = binarizeImagePair(nii1.img, nii2.img, threshold);
end

% Calculate overlap using image1 as reference
overlap1 = sum(nii1.img(:) .* nii2.img(:));
normOverlap1 = overlap1/sum(nii1.img(:));
if isnan(normOverlap1)
    normOverlap1 = 0;
end

% Reslice image1 using image2 as reference
reslicedImage1 = strrep(image1, '.nii', '_resliced.nii');
ea_fsl_reslice(image1, image2, reslicedImage1, 'trilinear', 0);

nii1 = ea_load_nii(reslicedImage1);
nii2 = ea_load_nii(image2);

% Binarize when needed
if binary
    [nii1.img, nii2.img] = binarizeImagePair(nii1.img, nii2.img, threshold);
end

% Calculate overlap using image2 as reference
overlap2 = sum(nii1.img(:) .* nii2.img(:));
normOverlap2 = overlap2/sum(nii2.img(:));
if isnan(normOverlap2)
    normOverlap2 = 0;
end

% Cleanup
ea_delete({reslicedImage1, reslicedImage2});


function [img1, img2] = binarizeImagePair(img1, img2, threshold)
% Binarize the image pair if they are not yet binaried

if (numel(unique(img1(:)))~=2 || numel(unique(img2(:)))~=2) ...
   || (all(unique(img1(:))==[0 1]') || all(unique(img2(:))==[0 1]'))
    if isnan(threshold)
        error('The images are not all binarized but threshold parameter is not supplied!');
    else
        if numel(threshold) == 1
            img1 = img1>threshold;
            img2 = img2>threshold;
        elseif numel(threshold) == 2
            img1 = img1>threshold(1);
            img2 = img2>threshold(2);
        end
    end
end

img1 = double(img1);
img2 = double(img2);
