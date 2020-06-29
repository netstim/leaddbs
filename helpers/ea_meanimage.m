function ea_meanimage(input, output)
% Calculate mean of input volumes
%
% Output directory is the same as the first input volume.
% Default output name is based on the name of the first input volume.
%
% Example:
%     ea_meanimage('A.nii') -> 'meanA.nii'
%     ea_meanimage({'A.nii,1', 'A.nii,2'}, 'mean') -> 'mean.nii'
%     ea_meanimage({'A.nii', 'B.nii'}) -> 'meanA.nii'
%     ea_meanimage({'A.nii,1', 'B.nii,1'}, 'mean') -> 'mean.nii'
%     ea_meanimage({'A.nii', 'B.nii,2'}, 'mean') -> 'mean.nii'

if ischar(input)
    input = ea_appendVolNum(input, 'all');
elseif iscell(input)
    if numel(input) ~= size(input, 1)
        input = input';
    end
    input = cellfun(@(X) ea_appendVolNum(X), input, 'Uni', 0);
end

[niifpath, niiname] = ea_niifileparts(input{1});

if nargin < 2
    output = ['mean', niiname];
end

outdir = fileparts(niifpath);

matlabbatch{1}.spm.util.imcalc.input = input;
matlabbatch{1}.spm.util.imcalc.output = output;
matlabbatch{1}.spm.util.imcalc.outdir = {outdir};
matlabbatch{1}.spm.util.imcalc.expression = 'mean(X)';
matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
matlabbatch{1}.spm.util.imcalc.options.dmtx = 1;
matlabbatch{1}.spm.util.imcalc.options.mask = 0;
matlabbatch{1}.spm.util.imcalc.options.interp = 1;
matlabbatch{1}.spm.util.imcalc.options.dtype = 4;

spm_jobman('run', {matlabbatch});
clear matlabbatch
