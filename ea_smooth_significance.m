function niicspmsig=ea_smooth_significance(XYZV,PTb,Vol,niic,options)
% function that generates a receptive field from datapoints by smoothing
% them. Then calculates an SPM statistic to check for significance.


odir=[options.root,options.patientname,filesep,'isosignificance',filesep];
mkdir(odir);
sdir=[odir,'SPM'];
try
    rmdir(sdir,'s');
end
mkdir(sdir)

XYZmm=[Vol.mat*[XYZV(:,1:3),ones(size(XYZV,1),1)]']';

% create bounding box version of files..
bb=[min(XYZmm(:,1)-5),min(XYZmm(:,2)-5),min(XYZmm(:,3)-5);
    max(XYZmm(:,1)+5),max(XYZmm(:,2)+5),max(XYZmm(:,3)+5)];
Vol.fname=[odir,'empty','.nii'];
spm_write_vol(Vol,niic);
ea_reslice_nii([odir,'empty','.nii'],[odir,'rempty','.nii'],[1,1,1]);
ea_crop_nii_bb([odir,'rempty','.nii'], bb);
%delete([odir,'empty','.nii']);

V=spm_vol([odir,'wrempty','.nii']);
niic=spm_read_vols(V);
XYZvx=V.mat\XYZmm';

XYZV(:,1:3)=round(XYZvx(1:3,:)');

XYZV(:,4)=zscore(XYZV(:,4));

for pt=1:max(PTb)
   niic(:)=nan;
    niic(sub2ind(size(niic),XYZV(PTb==pt,1),XYZV(PTb==pt,2),XYZV(PTb==pt,3)))=XYZV(PTb==pt,4);
    V.fname=[odir,num2str(pt),'.nii'];
    spm_write_vol(V,niic);
    spm_smooth(V.fname,V.fname,[3,3,3]);
    fis{pt}=[V.fname,',1'];
end

matlabbatch{1}.spm.stats.factorial_design.dir = {sdir};
matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = fis';
matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
spm_jobman('run',{matlabbatch});
clear matlabbatch


matlabbatch{1}.spm.stats.fmri_est.spmmat = {[sdir,filesep,'SPM.mat']};
matlabbatch{1}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;
spm_jobman('run',{matlabbatch});
clear matlabbatch

matlabbatch{1}.spm.stats.con.spmmat = {[sdir,filesep,'SPM.mat']};
matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'mainfx';
matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.delete = 1;
spm_jobman('run',{matlabbatch});
clear matlabbatch

matlabbatch{1}.spm.stats.results.spmmat = {[sdir,filesep,'SPM.mat']};
matlabbatch{1}.spm.stats.results.conspec.titlestr = '';
matlabbatch{1}.spm.stats.results.conspec.contrasts = 1;
matlabbatch{1}.spm.stats.results.conspec.threshdesc = 'FWE';
matlabbatch{1}.spm.stats.results.conspec.thresh = 0.05;
matlabbatch{1}.spm.stats.results.conspec.extent = 0;
matlabbatch{1}.spm.stats.results.conspec.mask.none = 1;
matlabbatch{1}.spm.stats.results.units = 1;
matlabbatch{1}.spm.stats.results.print = false;
matlabbatch{1}.spm.stats.results.write.tspm.basename = 'result';
spm_jobman('run',{matlabbatch});
clear matlabbatch

matlabbatch{1}.spm.util.imcalc.input = {Vol.fname;
    [sdir,filesep,'spmT_0001_result.nii,1']};
matlabbatch{1}.spm.util.imcalc.output = 'result.nii';
matlabbatch{1}.spm.util.imcalc.outdir = {''};
matlabbatch{1}.spm.util.imcalc.expression = 'i2';
matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
matlabbatch{1}.spm.util.imcalc.options.mask = 0;
matlabbatch{1}.spm.util.imcalc.options.interp = 1;
matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
spm_jobman('run',{matlabbatch});
clear matlabbatch




Vol2=spm_vol('result.nii');
niicspmsig=spm_read_vols(Vol);

for pt=1:max(PTb)
delete([odir,num2str(pt),'.nii']);
end
pth=fileparts(sdir);
rmdir(pth,'s');
