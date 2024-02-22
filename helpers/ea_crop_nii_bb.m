function ea_crop_nii_bb(input, bbox, output, interp)
% Crops nifti to a specific bounding-box.
arguments
    input   {mustBeFile}
    bbox    {mustBeNumeric}
    output  {mustBeTextScalar} = ''
    interp  {mustBeNumeric, mustBeMember(interp, [0, 1])} = 0 % 0 for nearest neighbor, 1 for trilinear
end

% Validate bbox size
if ~isequal(size(bbox), [2 3])
    ea_cprintf('CmdWinErrors', 'Wrong bounding box size! Should be 2x3.\n');
    return;
end

addpath(fullfile(fileparts(which('spm')), 'toolbox', 'OldNorm'));

if endsWith(input, '.gz')
    gzInput = 1;
    gunzip(input);
    input = erase(input, ".gz" + textBoundary('end'));
else
    gzInput = 0;
end

if endsWith(output, {'.nii', '.nii.gz'}) % output is file path
    prefixMode = 0;
    prefix = [ea_genid_rand, '_'];
    [inputFolder, inputName] = fileparts(input);
    tmpOutput = fullfile(inputFolder, [prefix, inputName, '.nii']);
    ea_mkdir(fileparts(GetFullPath(output)));
elseif isempty(output) % Overwrite input file
    prefixMode = 1;
    prefix = '';
    output = input;
else % output is prefix
    prefixMode = 1;
    prefix = output;
    [outputFolder, outputName] = fileparts(input);
    output = fullfile(outputFolder, [prefix, outputName, '.nii']);
end

V = spm_vol(input);

P = spm_imatrix(V(1).mat);
vox = P(7:9);
if any(vox<0)
    ea_reslice_nii(input, input, abs(vox), 0);
    V = spm_vol(input);
    P = spm_imatrix(V(1).mat);
    vox = P(7:9);
end

% create sn structure:
sn.VG = V;
sn.VF = V;
sn.Affine = eye(4);
sn.Tr = [];
sn.flags = [];

% create ropts structure:
ropts.preserve = 0;
ropts.bb = bbox;
ropts.vox = vox;
ropts.interp = interp;
ropts.wrap = [0,0,0];
ropts.prefix = prefix;
spm_write_sn([input,',1'], sn, ropts);

if prefixMode
    if gzInput
        gzip(output);
        delete(output);
    end
else
    if endsWith(output, '.nii') && ~strcmp(tmpOutput, output)
        movefile(tmpOutput, output);
    elseif endsWith(output, '.gz')
        gzip(tmpOutput);
        delete(tmpOutput);
        movefile([tmpOutput, '.gz'], output);
    end
end

% Delete unzipped input
if gzInput
    ea_delete(input);
end
