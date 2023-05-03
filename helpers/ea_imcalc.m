function nii = ea_imcalc(input, reference, prefix, interp, expr, maskFlag)
% warpper of SPM ImCalc

if ~exist('reference', 'var') || isempty(reference)
    reference = [ea_space, 't1.nii'];
end

% 'prefix' can be empty to overwrite the input image,
% or a prefix to prefix to input image (like in SPM, default is 'r'),
% or the real path of the desired output image.
if ~exist('prefix', 'var')
    prefix = 'r';
elseif endsWith(prefix, {'.nii', '.nii.gz'}) % 'prefix' is the real path of the output image
    outputImage = prefix;
    prefix = 'r';
end

if ~exist('interp', 'var')
    interp = 1; % use trilinear interpolation by default
end

if ~exist('expr', 'var')
    expr = 'i2';
end

if contains(expr, 'X')  % X
    dmtxFlag = 1;
else   % i1, i2, i3, ...
    dmtxFlag = 0;
end

if ~exist('maskFlag', 'var')
    % 0:  non implicit zero mask
    % 1:  implicit zero mask
    % -1: NaNs should be zeroed
    maskFlag = 0;
end

if strcmp(input(end-2:end),'.gz')
    gzinput = 1;
    gunzip(input);
    input = input(1:end-3);
else
    gzinput = 0;
end

input = GetFullPath(input);

if strcmp(reference(end-2:end),'.gz')
    gzreference = 1;
    gunzip(reference);
    reference = reference(1:end-3);
else
    gzreference = 0;
end

reference = GetFullPath(reference);

[fpath, fname]=ea_niifileparts(input);

matlabbatch{1}.spm.util.imcalc.input = {
                                        [reference,',1']
                                        [input,',1']
                                        };
matlabbatch{1}.spm.util.imcalc.output = [prefix, fname];
matlabbatch{1}.spm.util.imcalc.outdir = {fileparts(fpath)};
matlabbatch{1}.spm.util.imcalc.expression = expr;
matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
matlabbatch{1}.spm.util.imcalc.options.dmtx = dmtxFlag;
matlabbatch{1}.spm.util.imcalc.options.mask = maskFlag;
matlabbatch{1}.spm.util.imcalc.options.interp = interp;
matlabbatch{1}.spm.util.imcalc.options.dtype = 4;

spm_jobman('run',{matlabbatch});
clear matlabbatch

if nargout == 1
    nii = ea_load_nii([fileparts(fpath), filesep, prefix, fname, '.nii']);
end

if exist('outputImage', 'var')
    if endsWith(outputImage, '.nii')
        movefile([fileparts(fpath), filesep, prefix, fname, '.nii'], outputImage);
    elseif endsWith(outputImage, '.nii.gz')
        outputImage = regexprep(outputImage, '\.nii\.gz$', '\.nii');
        movefile([fileparts(fpath), filesep, prefix, fname, '.nii'], outputImage);
        gzip(outputImage);
        delete(outputImage);
    end
else
    outputImage = [fileparts(fpath), filesep, prefix, fname, '.nii'];
    if gzinput
        gzip(outputImage);
        delete(outputImage);
    end
end

if gzinput
    gzip(input);
    delete(input);
end

if gzreference
    delete(reference);
end
