nii=load_nii('Thal.nii')
[xx,yy,zz]=ind2sub(size(nii.img),find(nii.img));
copt.cut_from_L=size(nii.img,1)-(max(xx));
copt.cut_from_R=(min(xx));
copt.cut_from_A=size(nii.img,2)-(max(yy));
copt.cut_from_P=min(yy);
copt.cut_from_I=size(nii.img,3)-(max(zz));
copt.cut_from_S=min(zz);

nii=clip_nii(nii,copt);
save_nii(nii,'cThal.nii')
