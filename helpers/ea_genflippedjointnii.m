function joints=ea_genflippedjointnii(right, left,keepbinary)
% function to create a joint flipped nifti roi e.g. based on left and right
% vtas.

if ~exist('keepbinary','var')
    keepbinary=0;
end
[rpth,rfn,rext]=fileparts(right);
[lpth,lfn,lext]=fileparts(left);
joints={fullfile(rpth,['bihem_',rfn,rext]),fullfile(lpth,['bihem_',lfn,lext])};

if (~exist(fullfile(rpth,['bihem_',rfn,rext]),'file')) || (~exist(fullfile(lpth,['bihem_',lfn,lext]),'file')) 
    ea_flip_lr_nonlinear(right,fullfile(rpth,['fl_',rfn,rext]),0);
    ea_flip_lr_nonlinear(left,fullfile(lpth,['fl_',lfn,lext]),0);
    jright=ea_load_nii(right);
    flleft=ea_load_nii(fullfile(lpth,['fl_',lfn,lext]));
    jright.img=(jright.img+flleft.img)/2;
    jright.fname=fullfile(rpth,['bihem_',rfn,rext]);
    if keepbinary
        jright.img=jright.img>0;
    end
    ea_write_nii(jright);
    
    jleft=ea_load_nii(left);
    flright=ea_load_nii(fullfile(rpth,['fl_',rfn,rext]));
    jleft.img=(jleft.img+flright.img)/2;
    jleft.fname=fullfile(lpth,['bihem_',lfn,lext]);
    if keepbinary
        jleft.img=jleft.img>0;
    end
    ea_write_nii(jleft);
end




