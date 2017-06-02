function ea_autocrop(image, prefix, mask)
% Crop the image to its minimum bounding box

if nargin < 2
    prefix = '';	% overwrite the image by default
end

if nargin < 3
    mask = 0;	% do not mask the image (remove background) by default
end

[bbox, BW] = ea_autobbox(image);

if mask
    nii = load_nii(image);
    nii.img = nii.img .* BW;

    if isempty(prefix)
        output = image;
    else
        [pth, fname, ext] = fileparts(image);
        output = fullfile(pth,[prefix, fname, ext]);
    end

    save_nii(nii, output);
    ea_crop_nii_bb(output, '', bbox);
else
    ea_crop_nii_bb(image, prefix, bbox);
end
