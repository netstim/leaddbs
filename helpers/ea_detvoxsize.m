function vsize=ea_detvoxsize(mat)
% Determine the voxel size based on the affine matrix of the image
% 'mat' can be either the affine matrix or the nifti file path

vox = [1,1,1;
       1,1,1;
       1,1,1;
       2,1,1;
       1,2,1;
       1,1,2];
mm = ea_vox2mm(vox, mat);
try
    vsize = diag(ea_pdist2(mm(1:3,:),mm(4:6,:)))';
catch % auto load lead path.
    lead path
    vsize = diag(ea_pdist2(mm(1:3,:),mm(4:6,:)))';
end
