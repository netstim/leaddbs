function ea_split_nii_lr(fname,offset)
% offset may be used to not split at x = 0 mm but at abs(x)>offset.

if ~exist('offset','var')
    offset=0;
end

nii=ea_load_nii(fname);
lnii=nii;
rnii=nii;

if ~isBIDSFileName(fname)
    [pth,fn,ext]=fileparts(fname);
    lnii.fname=fullfile(pth,[fn,'_l',ext]);
    rnii.fname=fullfile(pth,[fn,'_r',ext]);
else
    lnii.fname = setBIDSEntity(fname, 'hemi', 'L');
    rnii.fname = setBIDSEntity(fname, 'hemi', 'R');
end

[xx,yy,zz]=ind2sub(size(nii.img),1:numel(nii.img));
XYZ=[xx;yy;zz;ones(1,length(xx))];
XYZ=nii.mat*XYZ;

lnix=XYZ(1,:)<offset;
rnii.img(lnix)=0;

rnix=XYZ(1,:)>-offset;
lnii.img(rnix)=0;

ea_write_nii(rnii);
ea_write_nii(lnii);
