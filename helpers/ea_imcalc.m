function nii = ea_imcalc(input, output, expr, opts)
% warpper of SPM ImCalc

arguments
    input   {mustBeText}
    output  {mustBeTextScalar}
    expr    {mustBeTextScalar} = 'i2' % Reslice the image to to match reference by default
    opts.mask    {mustBeNumeric} = 0 % 0: non implicit zero mask, 1:  implicit zero mask, -1: NaNs should be zeroed
    opts.interp  {mustBeNumeric} = 1  % Use trilinear interpolation by default
    opts.dtype   {mustBeNumeric} = 4  % 2: 'uint8', 4: 'int16', 8: 'int32', 16: 'float32', 64: 'float64', 256: 'int8', 512: 'uint16', 768: 'uint32'
    opts.matchBBox {mustBeNumericOrLogical} = false % Do not match the bbox of inputs by default
end

% Make sure input to SPM ImCalc is cell
if ischar(input)
    input = {input};
end

% Make sure input to SPM ImCalc is column vector
if isrow(input)
    input = input';
end

% Get full path of the input
input = GetFullPath(input);

% Unzip input incase necessray
gzInputs = input(endsWith(input, '.gz'));
if ~isempty(gzInputs)
    gunzip(gzInputs);
    input = erase(input, ".gz" + textBoundary('end'));
end

% Match bounding box
if opts.matchBBox
    % Calculate largest bounding box
    [bbox, voxsize] = cellfun(@(x) spm_get_bbox(x, 'nz'), input, 'Uni', 0);
    bbox = cell2mat(bbox);
    bbox = [min(bbox(1:2:end, :)); max(bbox(2:2:end, :))];

    % Crop input to the largest bounding boxs
    tempFolder = [ea_getleadtempdir, ea_genid_rand(1,6)];
    ea_mkdir(tempFolder);
    [~, inputNames] = fileparts(input);
    croppedFile = strcat(tempFolder, filesep, inputNames, '.nii');
    cellfun(@(x, y) ea_crop_nii_bb(x, bbox, y), input, croppedFile);

    % Reslice the first image to the finest resolution
    ea_resample_image_by_spacing(croppedFile{1}, min(abs(cell2mat(voxsize))), 0, 0, ~opts.interp, croppedFile{1});

    % Override input
    input = croppedFile;
end

% Parse output
output = GetFullPath(output);
[outputWithoutExt, outputFileName] = ea_niifileparts(output);

if contains(expr, 'X')  % X
    dmtxFlag = 1;
else   % i1, i2, i3, ...
    dmtxFlag = 0;
end

% Run SPM ImCalc
matlabbatch{1}.spm.util.imcalc.input = strcat(input, ',1');
matlabbatch{1}.spm.util.imcalc.output = outputFileName;
matlabbatch{1}.spm.util.imcalc.outdir = {fileparts(outputWithoutExt)};
matlabbatch{1}.spm.util.imcalc.expression = expr;
matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
matlabbatch{1}.spm.util.imcalc.options.dmtx = dmtxFlag;
matlabbatch{1}.spm.util.imcalc.options.mask = opts.mask;
matlabbatch{1}.spm.util.imcalc.options.interp = opts.interp;
matlabbatch{1}.spm.util.imcalc.options.dtype = opts.dtype;

spm_jobman('run', {matlabbatch});
clear matlabbatch

if nargout == 1
    nii = ea_load_nii([outputWithoutExt, '.nii']);
end

% Zip in case .gz output specified
if endsWith(output, '.gz')
    gzip([outputWithoutExt, '.nii']);
    delete([outputWithoutExt, '.nii']);
end

% Delete unzipped input
if ~isempty(gzInputs)
    ea_delete(erase(gzInputs, ".gz" + textBoundary('end')));
end

% Delete cropped image
if opts.matchBBox
    ea_delete(tempFolder);
end
