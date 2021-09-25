function ea_genflippedjointnii(right,left)
% function to flip nifti roi e.g. based on left and right
% vtas.

[rpth,rfn,rext]=fileparts(right);
[lpth,lfn,lext]=fileparts(left);
if ~exist(fullfile(rpth,['fl_',rfn,rext]),'file')
    try
        ea_flip_lr_nonlinear(right,fullfile(rpth,['fl_',rfn,rext]),0);
    catch
        ea_reslice_nii(right,fullfile(rpth,['fl_',rfn,rext]));
        ea_flip_lr_nonlinear(fullfile(rpth,['fl_',rfn,rext]),fullfile(rpth,['fl_',rfn,rext]),0);
    end
end
if ~exist(fullfile(lpth,['fl_',lfn,lext]),'file')
    try
        ea_flip_lr_nonlinear(left,fullfile(lpth,['fl_',lfn,lext]),0);
    catch
        ea_reslice_nii(left,fullfile(rpth,['fl_',lfn,rext]));
        ea_flip_lr_nonlinear(fullfile(rpth,['fl_',lfn,rext]),fullfile(rpth,['fl_',lfn,rext]),0);     
    end
end