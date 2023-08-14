% fixes OASIS TPM for use in SPM - based on post by J. A. https://www.jiscmail.ac.uk/cgi-bin/webadmin?A2=spm;a6060c5a.1609
% THIS HAS ALREADY BEEN APPLIED TO THE FILE PRESENT IN YOUR FOLDER AND IS
% ONLY MEANT FOR REPRODUCIBILITY REASONS. PLEASE DON'T RUN THIS FILE.
copyfile('TPM.nii','TPM_backup.nii');

matlabbatch{1}.spm.spatial.smooth.data = {'TPM.nii'};
matlabbatch{1}.spm.spatial.smooth.fwhm = [3 3 3];
matlabbatch{1}.spm.spatial.smooth.dtype = 0;
matlabbatch{1}.spm.spatial.smooth.im = 0;
matlabbatch{1}.spm.spatial.smooth.prefix = 's';
spm_jobman('run',{matlabbatch});
clear matlabbatch

nii=ea_load_untouch_nii('sTPM.nii');
nii.img=nii.img+(1/1000); % add constant
for c=1:size(nii.img,4) % read into nii.img matrix X
    thisc=nii.img(:,:,:,c);
   X(c,:)=thisc(:); 
end

X=X./repmat(sum(X,1),6,1); % make sure sum eq to 1 for each voxel
X=X*32767;
sum(X,1);

for c=1:size(nii.img,4) % write back into images
    thisc(:)=X(c,:);
   nii.img(:,:,:,c)=thisc; 
end


ea_save_untouch_nii(nii,'TPM.nii');

