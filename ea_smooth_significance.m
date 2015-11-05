function niicspmsig=ea_smooth_significance(XYZV,Vol,niic,options)
% function that generates a receptive field from datapoints by smoothing
% them. Then calculates an SPM statistic to check for significance.


odir=[options.root,options.patientname,filesep,'isosignificance',filesep];
mkdir(odir);
sdir=[odir,'SPM'];
mkdir(sdir)

XYZmm=[Vol.mat*[XYZV(:,1:3),ones(size(XYZV,1),1)]']';

% create bounding box version of files..
bb=[min(XYZmm(:,1)-5),min(XYZmm(:,2)-5),min(XYZmm(:,3)-5);
    max(XYZmm(:,1)+5),max(XYZmm(:,2)+5),max(XYZmm(:,3)+5)];
Vol.fname=[odir,'empty','.nii'];
spm_write_vol(Vol,niic);
ea_crop_nii_bb(Vol.fname,'w',bb);
delete([odir,'empty','.nii']);

V=spm_vol([odir,'wempty','.nii']);
XYZvx=V.mat\XYZmm'; 

XYZV(:,1:3)=XYZvx(1:3,:)';

XYZV(:,4)=zscore(XYZV(:,4));

for el=1:size(XYZV,1)
   niic(:)=nan;
   
    
    niic(XYZV(el,1),XYZV(el,2),XYZV(el,3))=XYZV(el,4);
    V.fname=[odir,num2str(el),'.nii'];
    spm_write_vol(V,niic);
    spm_smooth(V.fname,V.fname,[3,3,3]);
    fis{el}=[V.fname,',1'];
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
cfg_util('run',{matlabbatch});
clear matlabbatch


matlabbatch{1}.spm.stats.fmri_est.spmmat = {[sdir,filesep,'SPM.mat']};
matlabbatch{1}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;
cfg_util('run',{matlabbatch});
clear matlabbatch

matlabbatch{1}.spm.stats.con.spmmat = {[sdir,filesep,'SPM.mat']};
matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'mainfx';
matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.delete = 1;
cfg_util('run',{matlabbatch});
clear matlabbatch

keyboard