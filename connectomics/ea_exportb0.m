function ea_exportb0(options)

bvals=load([options.root,options.patientname,filesep,options.prefs.bval]);
idx=find(bvals<10);
cnt=1;

if size(idx,1)<size(idx,2)
    idx=idx';
end
for fi=idx'

   fis{cnt}=[options.root,options.patientname,filesep,options.prefs.dti,',',num2str(fi)];
   %nii=ea_load_nii(fis{cnt});
   %X(:,:,:,cnt)=nii.img;
   cnt=cnt+1;
end



if length(fis)==1
 expr='X';   
else
    expr='mean(X)';
end
    
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
