function ea_autocrop(input, output, opts)
% Crop the image to its minimum bounding box
arguments
    input       {mustBeFile}
    output      {mustBeTextScalar} = '' % Overwrite input by default
    opts.mask   {mustBeNumericOrLogical} = false % Do not mask the image (remove background) by default
    opts.margin {mustBeNumeric} = 5 % Add margin to the cropped image, unit in voxel
    opts.method {mustBeMember(opts.method, {'spm', 'multithresh'})} = 'spm' % Method to get the bbox
end

% Unzip input incase necessray
if endsWith(input, '.gz')
    wasGzip = 1;
    gunzip(input);
    input = erase(input, ".gz" + textBoundary('end'));
else
    wasGzip = 0;
end

switch opts.method
    case 'spm' % Use SPM to get the bbox: simple non-zero and non-nan thresholding
        bbox = spm_get_bbox(input, 'nz');
        if opts.margin
            affine = ea_get_affine(input);
            voxbbox = ea_mm2vox(bbox, affine);
            voxbbox(1,:) = voxbbox(1,:) - opts.margin;
            voxbbox(2,:) = voxbbox(2,:) + opts.margin;
            bbox = ea_vox2mm(voxbbox, affine);
        end
    case 'multithresh' % Use multithresh for thresholding
        [bbox, BW] = ea_autobbox(input, opts.margin);
end

if opts.mask
    nii = load_nii(input);
    nii.img(isnan(nii.img)) = 0; 
    if strcmp(opts.method, 'multithresh')
        nii.img = nii.img .* BW;
    end

    [inputFolder, inputName, inputExt] = fileparts(input);
    maskedInput = fullfile(inputFolder, [ea_genid_rand, '_', inputName, inputExt]);

    save_nii(nii, maskedInput);

    ea_crop_nii_bb(maskedInput, bbox, output);

    ea_delete(maskedInput);
else
    ea_crop_nii_bb(input, bbox, output);
end

if wasGzip
    delete(input);
end
