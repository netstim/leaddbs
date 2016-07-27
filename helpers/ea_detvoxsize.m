function vsize=ea_detvoxsize(mat)
% determine the voxel size

vox = [1,1,1;
       1,1,1;
       1,1,1;
       2,1,1;
       1,2,1;
       1,1,2];
mm = vox2mm(vox, mat);
vsize = diag(pdist2(mm(1:3,:),mm(4:6,:)))';
