function [overlap1, normOverlap1, overlap2, normOverlap2] = ea_nii_overlap(image1, image2, binary, threshold)
% Calculate the overlap between two images (within the same space)
%
% Return the overlaps in VOXEL (NOT IN MM^3), and normalized overlaps (
% normalized by the voxels of nii1 and nii2 respectively)

overlap1 = 0;
normOverlap1 = 0;
overlap2 = 0;
normOverlap2 = 0;

% Use binary calculation default
if ~exist('binary', 'var')
    binary = 1;
end

nii1 = load_untouch_nii(image1);
nii2 = load_untouch_nii(image2);

if binary
    % Binarize the image if they are not yet binaried
    if (numel(unique(nii1.img(:)))~=2 || numel(unique(nii2.img(:)))~=2) ...
       || (all(unique(nii1.img(:))==[0 1]') || all(unique(nii2.img(:))==[0 1]'))
        if ~exist('threshold', 'var')
            error('The images are not all binarized but threshold parameter is not supplied!');
        else
            if numel(threshold) == 1
                nii1.img = nii1.img>=threshold;
                nii2.img = nii2.img>=threshold;
            elseif numel(threshold) == 2
                nii1.img = nii1.img>=threshold(1);
                nii2.img = nii2.img>=threshold(2);
            end
        end
    end
end

nii1.img = double(nii1.img);
nii2.img = double(nii2.img);

% Map to non-zeros voxel in image2 to image1
[xvox, yvox, zvox] = ind2sub(size(nii2.img), find(nii2.img(:)));
mm = ea_vox2mm([xvox, yvox, zvox], image2);
vox = round(ea_mm2vox(mm, image1));
filter = all(vox>[0,0,0],2) & all(vox<=size(nii1.img),2); % Remove voxel out of bbox
if sum(filter)
    vox = vox(filter, :);
    ind = unique(sub2ind(size(nii1.img), vox(:,1), vox(:,2), vox(:,3)));
    % Checking overlap for image1
    overlap1 = sum(nii1.img(ind));
    normOverlap1 = overlap1/sum(nii1.img(:));
end

% Map to non-zeros voxel in image1 to image2
[xvox, yvox, zvox] = ind2sub(size(nii1.img), find(nii1.img(:)));
mm = ea_vox2mm([xvox, yvox, zvox], image1);
vox = round(ea_mm2vox(mm, image2));
filter = all(vox>[0,0,0],2) & all(vox<=size(nii2.img),2); % Remove voxel out of bbox
if sum(filter)
    vox = vox(filter, :);
    ind = unique(sub2ind(size(nii2.img), vox(:,1), vox(:,2), vox(:,3)));
    % Checking overlap for image2
    overlap2 = sum(nii2.img(ind));
    normOverlap2 = overlap2/sum(nii2.img(:));
end
