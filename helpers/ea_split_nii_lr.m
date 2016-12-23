function ea_split_nii_lr(fname)

nii=ea_load_nii(fname);
lnii=nii;
rnii=nii;
[pth,fn,ext]=fileparts(fname);
lnii.fname=fullfile(pth,[fn,'_l',ext]);
rnii.fname=fullfile(pth,[fn,'_r',ext]);

[xx,yy,zz]=ind2sub(size(nii.img),1:numel(nii.img));
XYZ=[xx;yy;zz;ones(1,length(xx))];
XYZ=nii.mat*XYZ;

lnix=XYZ(1,:)<0;
rnii.img(lnix)=0;

rnix=XYZ(1,:)>0;
lnii.img(rnix)=0;

ea_write_nii(rnii);
ea_write_nii(lnii);
