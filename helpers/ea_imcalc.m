function nii = ea_imcalc(input, reference, prefix, interp, expr)
% warpper of SPM ImCalc

if ~exist('reference', 'var') || isempty(reference)
    reference = [ea_space, 't1.nii'];
end

if ~exist('prefix', 'var')
    prefix = 'r';
end

if ~exist('interp', 'var')
    interp = 1; % use trilinear interpolation by default
end

if ~exist('expr', 'var')
    expr = 'i2';
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
matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
matlabbatch{1}.spm.util.imcalc.options.mask = 0;
matlabbatch{1}.spm.util.imcalc.options.interp = interp;
matlabbatch{1}.spm.util.imcalc.options.dtype = 4;

spm_jobman('run',{matlabbatch});
clear matlabbatch

if nargout == 1
    nii = ea_load_nii([fileparts(fpath), filesep, prefix, fname, '.nii']);
end

if gzinput
    gzip(input);
    delete(input);
end

if gzreference
    delete(reference);
end
