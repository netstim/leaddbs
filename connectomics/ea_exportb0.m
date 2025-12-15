function ea_exportb0(options)

b0threshold = 10;

disp('Export b0...');
bvals=load([options.root,options.patientname,filesep,options.prefs.bval]);
idx=find(bvals<b0threshold);

if isempty(idx)
    matlab.desktop.editor.openAndGoToLine(which('ea_exportb0'), 3);
    error(sprintf(['Exporting b0 image failed: b0 threshold is too small!\n' ...
           'Please check your dti.bval and _temporarily_ set ''b0threshold''', ...
           ' in function ''ea_exportb0'' to a proper higher value.']));
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
    % BIDS: path includes subdirectory (e.g., 'dwi/sub-..._b0.nii')
    outdir = fileparts([options.root, options.patientname, filesep, options.prefs.b0]);
else
    % Classic: just filename
    outdir = [options.root, options.patientname];
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

% BIDS FIX: Correct b0 header origin if it's far from image center
% This fixes issues with DWI data that have incorrect qform/sform origins
b0Path = fullfile(outdir, [b0Name, '.nii']);
try
    V = spm_vol(b0Path);
    origin_vox = [-V.mat(1,4)/V.mat(1,1), -V.mat(2,4)/V.mat(6), -V.mat(3,4)/V.mat(11)];
    center_vox = V.dim / 2;
    
    % If origin is more than 20% off from center, correct it
    if abs(origin_vox(3) - center_vox(3)) > 0.2 * center_vox(3)
        fprintf('\nCorrecting b0 header origin (was at voxel [%.1f, %.1f, %.1f], center is [%.1f, %.1f, %.1f])...\n', ...
            origin_vox(1), origin_vox(2), origin_vox(3), center_vox(1), center_vox(2), center_vox(3));
        
        % Read the image
        img_data = spm_read_vols(V);
        
        % Create new header with corrected origin (centered)
        V_new = V;
        V_new.mat(1,4) = -V.mat(1,1) * center_vox(1);
        V_new.mat(2,4) = -V.mat(6) * center_vox(2);
        V_new.mat(3,4) = -V.mat(11) * center_vox(3);
        
        % Write with corrected header
        spm_write_vol(V_new, img_data);
        fprintf('b0 header corrected successfully.\n\n');
    end
catch ME
    warning('Failed to correct b0 header: %s', ME.message);
end

