function ea_crop_nii_bb(input, bbox, output)
% Crops nifti to a specific bounding-box.
arguments
    input   {mustBeFile}
    bbox    {mustBeNumeric} % in mm
    output  {mustBeTextScalar} = ''
end

% Validate bbox size
if ~isequal(size(bbox), [2 3])
    ea_cprintf('CmdWinErrors', 'Wrong bounding box size! Should be 2x3.\n');
    return;
end

% Check output parameter
input = GetFullPath(input);
if endsWith(output, {'.nii', '.nii.gz'}) % output is file path
    output = GetFullPath(output);
    ea_mkdir(fileparts(output));
elseif isempty(output) % Overwrite input file
    output = input;
else % prefix input (output is prefix)
    [inputFolder, inputName, inputExt] = fileparts(input);
    output = fullfile(inputFolder, [output, inputName, inputExt]);
end

% Determine output format
if endsWith(output, '.nii')
    envStr = 'FSLOUTPUTTYPE=NIFTI';
elseif endsWith(output, '.nii.gz')
    envStr = 'FSLOUTPUTTYPE=NIFTI_GZ';
end

% Get voxel bbox to be used for fslroi (indexing starts with 0)
vbbox = round(ea_mm2vox(bbox, input, 0));

% Construct [xyz]min and [xyz]size to be used for fslroi
xmin = num2str(vbbox(1, 1));
xsize = num2str(vbbox(2, 1) - vbbox(1, 1) + 1);
ymin = num2str(vbbox(1, 2));
ysize = num2str(vbbox(2, 2) - vbbox(1, 2) + 1);
zmin = num2str(vbbox(1, 3));
zsize = num2str(vbbox(2, 3) - vbbox(1, 3) + 1);

% Get fslroi binary path
fslroi = ea_getExec(fullfile(ea_getearoot, 'ext_libs', 'fsl', 'fslroi'), escapePath=true);

% Construct cmd
cmd = [fslroi ' ' ea_path_helper(input) ' ' ea_path_helper(output) ' ' xmin ' ' xsize ' ' ymin ' ' ysize ' ' zmin ' ' zsize];

% Run fslroi to crop image
ea_runcmd(cmd, env=envStr);
