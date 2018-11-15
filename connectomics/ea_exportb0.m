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

if length(fis)==1
    expr='i1';

    matlabbatch{1}.spm.util.imcalc.input = fis';
    matlabbatch{1}.spm.util.imcalc.output = [options.prefs.b0];
    matlabbatch{1}.spm.util.imcalc.outdir = {[options.root,options.patientname]};
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
    matlabbatch{1}.spm.util.imcalc.output = [options.prefs.b0];
    matlabbatch{1}.spm.util.imcalc.outdir = {[options.root,options.patientname]};
    matlabbatch{1}.spm.util.imcalc.expression = expr;
    matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
    matlabbatch{1}.spm.util.imcalc.options.dmtx = 1;
    matlabbatch{1}.spm.util.imcalc.options.mask = 0;
    matlabbatch{1}.spm.util.imcalc.options.interp = 1;
    matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
    spm_jobman('run',{matlabbatch}); clear matlabbatch
end

