function ea_nc_segment(options)

%% neck-crop and segment MR-image
% check if pre_tra has already been segmented..

matlabbatch{1}.spm.tools.preproc8.channel.vols = {[options.root,options.patientname,filesep,options.prefs.prenii_unnormalized,',1']};
matlabbatch{1}.spm.tools.preproc8.channel.biasreg = 0.0001;
matlabbatch{1}.spm.tools.preproc8.channel.biasfwhm = 60;
matlabbatch{1}.spm.tools.preproc8.channel.write = [0 0];
matlabbatch{1}.spm.tools.preproc8.tissue(1).tpm = {[options.earoot,'templates',filesep,'TPM.nii,1']};
matlabbatch{1}.spm.tools.preproc8.tissue(1).ngaus = 2;
matlabbatch{1}.spm.tools.preproc8.tissue(1).native = [1 0];
matlabbatch{1}.spm.tools.preproc8.tissue(1).warped = [0 0];
matlabbatch{1}.spm.tools.preproc8.tissue(2).tpm = {[options.earoot,'templates',filesep,'TPM.nii,2']};
matlabbatch{1}.spm.tools.preproc8.tissue(2).ngaus = 2;
matlabbatch{1}.spm.tools.preproc8.tissue(2).native = [1 0];
matlabbatch{1}.spm.tools.preproc8.tissue(2).warped = [0 0];
matlabbatch{1}.spm.tools.preproc8.tissue(3).tpm = {[options.earoot,'templates',filesep,'TPM.nii,3']};
matlabbatch{1}.spm.tools.preproc8.tissue(3).ngaus = 2;
matlabbatch{1}.spm.tools.preproc8.tissue(3).native = [1 0];
matlabbatch{1}.spm.tools.preproc8.tissue(3).warped = [0 0];
matlabbatch{1}.spm.tools.preproc8.tissue(4).tpm = {[options.earoot,'templates',filesep,'TPM.nii,4']};
matlabbatch{1}.spm.tools.preproc8.tissue(4).ngaus = 3;
matlabbatch{1}.spm.tools.preproc8.tissue(4).native = [1 0];
matlabbatch{1}.spm.tools.preproc8.tissue(4).warped = [0 0];
matlabbatch{1}.spm.tools.preproc8.tissue(5).tpm = {[options.earoot,'templates',filesep,'TPM.nii,5']};
matlabbatch{1}.spm.tools.preproc8.tissue(5).ngaus = 4;
matlabbatch{1}.spm.tools.preproc8.tissue(5).native = [1 0];
matlabbatch{1}.spm.tools.preproc8.tissue(5).warped = [0 0];
matlabbatch{1}.spm.tools.preproc8.tissue(6).tpm = {[options.earoot,'templates',filesep,'TPM.nii,6']};
matlabbatch{1}.spm.tools.preproc8.tissue(6).ngaus = 2;
matlabbatch{1}.spm.tools.preproc8.tissue(6).native = [0 0];
matlabbatch{1}.spm.tools.preproc8.tissue(6).warped = [0 0];
matlabbatch{1}.spm.tools.preproc8.warp.mrf = 0;
matlabbatch{1}.spm.tools.preproc8.warp.reg = 4;
matlabbatch{1}.spm.tools.preproc8.warp.affreg = 'mni';
matlabbatch{1}.spm.tools.preproc8.warp.samp = 3;
matlabbatch{1}.spm.tools.preproc8.warp.write = [0 0];
spm_jobman('run',{matlabbatch});
clear matlabbatch



Vskull=spm_vol([options.root,options.patientname,filesep,'c4',options.prefs.prenii_unnormalized]);

Xskull=spm_read_vols(Vskull);
[xx,yy,zz]=ind2sub(size(Xskull),find(Xskull>0.2));
    delete([options.root,options.patientname,filesep,'c4',options.prefs.prenii_unnormalized]);


% crop:
V=spm_vol([options.root,options.patientname,filesep,options.prefs.prenii_unnormalized]);
X=spm_read_vols(V);




% skinstrip:
M=Xskull>0.2;
clear Xskull

for c=1:3
    Vm=spm_vol([options.root,options.patientname,filesep,'c',num2str(c),options.prefs.prenii_unnormalized]);
    Xm=spm_read_vols(Vm);
    delete([options.root,options.patientname,filesep,'c',num2str(c),options.prefs.prenii_unnormalized]);
    M=M+Xm>0.1;

end
clear Xm

X=X.*M;

X(max(xx):end,:,:)=0;
X(:,max(yy):end,:)=0;
X(:,:,max(zz):end)=0;
X(1:min(xx),:,:)=0;
X(:,1:min(yy),:)=0;
X(:,:,1:min(zz))=0;

V.fname=[options.root,options.patientname,filesep,'c',options.prefs.prenii_unnormalized];
spm_write_vol(V,X);
