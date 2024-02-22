function ea_autocrop(input, output, mask, margin)
% Crop the image to its minimum bounding box
arguments
    input   {mustBeFile}
    output  {mustBeTextScalar} = '' % Overwrite input by default
    mask    {mustBeNumericOrLogical} = false % Do not mask the image (remove background) by default
    margin  {mustBeNumeric} = 5 % Add margin to the cropped image, unit in voxel
end

[bbox, BW] = ea_autobbox(input, margin);

if mask
    nii = load_nii(input);
    nii.img = nii.img .* BW;

    [inputFolder, inputName, inputExt] = fileparts(input);
    maskedInput = fullfile(inputFolder, [ea_genid_rand, '_', inputName, inputExt]);

    save_nii(nii, maskedInput);

    ea_crop_nii_bb(maskedInput, bbox, output);

    ea_delete(maskedInput);
else
    ea_crop_nii_bb(input, bbox, output);
end
