function ea_exportb0(options)

disp('Export b0...');
bvals=load([options.root,options.patientname,filesep,options.prefs.bval]);
b0threshold = ea_detect_b0threshold(bvals);
idx=find(bvals<=b0threshold);

if isempty(idx)
    error(['Exporting b0 image failed: no volumes detected below the auto b0 threshold (', ...
           num2str(b0threshold), '). Please check your dti.bval file.']);
end

cnt=1;

if size(idx,1)<size(idx,2)
    idx=idx';
end
for fi=idx'
   fis{cnt}=[options.root,options.patientname,filesep,options.prefs.dti,',',num2str(fi)];
   cnt=cnt+1;
end

% Determine output directory (handle BIDS paths with subdirectories)
[~, b0Name] = fileparts(options.prefs.b0);
if contains(options.prefs.b0, filesep)
    % BIDS: path includes subdirectory (e.g., 'preprocessing/dwi/sub-..._b0.nii')
    outdir = fileparts(fullfile(options.root, options.patientname, options.prefs.b0));
else
    % Classic: just filename
    outdir = fullfile(options.root, options.patientname);
end

if length(fis)==1
    expr='i1';

    matlabbatch{1}.spm.util.imcalc.input = fis';
    matlabbatch{1}.spm.util.imcalc.output = [b0Name, '.nii'];
    matlabbatch{1}.spm.util.imcalc.outdir = {outdir};
    matlabbatch{1}.spm.util.imcalc.expression = expr;
    matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
    matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
    matlabbatch{1}.spm.util.imcalc.options.mask = 0;
    matlabbatch{1}.spm.util.imcalc.options.interp = 1;
    matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
    spm_jobman('run',{matlabbatch}); clear matlabbatch
else
    expr='mean(X)';

    matlabbatch{1}.spm.util.imcalc.input = fis';
    matlabbatch{1}.spm.util.imcalc.output = [b0Name, '.nii'];
    matlabbatch{1}.spm.util.imcalc.outdir = {outdir};
    matlabbatch{1}.spm.util.imcalc.expression = expr;
    matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
    matlabbatch{1}.spm.util.imcalc.options.dmtx = 1;
    matlabbatch{1}.spm.util.imcalc.options.mask = 0;
    matlabbatch{1}.spm.util.imcalc.options.interp = 1;
    matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
    spm_jobman('run',{matlabbatch}); clear matlabbatch
end


function thr = ea_detect_b0threshold(bvals)
% Use the lowest b-value shell as b0 and separate it from the next shell.

bvals = bvals(:);
sv = sort(unique(bvals));

if numel(sv) < 2
    thr = sv(1) + 1;
    return;
end

% Minimum gap required to treat the next b-value as a new shell
% Smaller gaps are absorbed into the current low-b shell
shellTol = 100;

% Grow the low-b shell while neighboring b-values remain close
topIdx = 1;
while topIdx < numel(sv) && (sv(topIdx+1) - sv(topIdx)) <= shellTol
    topIdx = topIdx + 1;
end

if topIdx == numel(sv)
    % No higher shell to separate from
    thr = sv(end) + 1;
    return;
end

% Select everything up to the low shell, but not the next shell
thr = (sv(topIdx) + sv(topIdx+1)) / 2;
