function ea_largestcomponent_nii(fname,cnum)
% helper function that reduces the nifti components to the largest
% component.

nii=ea_load_nii(fname);
nii.img(isnan(nii.img))=0;
nii.img=logical(nii.img);
C=bwconncomp(nii.img);


ls=cellfun(@length,C.PixelIdxList,'Uniformoutput',0);
ls=cell2mat(ls);
[~,ix]=sort(ls);

[pth,fn,ext]=fileparts(fname);

for cc=0:cnum-1
   nii.img(:)=0;
   nii.img(C.PixelIdxList{ix(end-cc)})=1;
   nii.fname=fullfile(pth,[fn,'_c',num2str(cc),ext]);
   spm_write_vol(nii,nii.img);
   
end