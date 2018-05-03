function ea_symmetrize_nii_pair(p1,p2,interp)

if ~exist('interp','var')
    interp=1;
end
[pth,fn,ext]=fileparts(p1);
po1=fullfile(pth,[fn,'_flipped',ext]);
[pth,fn,ext]=fileparts(p2);
po2=fullfile(pth,[fn,'_flipped',ext]);
ea_flip_lr_nonlinear(p1,po1,interp);
ea_flip_lr_nonlinear(p2,po2,interp);

ori1=ea_load_nii(p1);
ori2=ea_load_nii(p2);

ea_conformspaceto(p1,po2,interp,[],[],0);
ea_conformspaceto(p2,po1,interp,[],[],0);
fl1=ea_load_nii(po1);
fl2=ea_load_nii(po2);

sym1=ori1;
sym1.img=(ori1.img+fl2.img)/2;
sym2=ori2;
sym2.img=(ori2.img+fl1.img)/2;

[pth,fn,ext]=fileparts(sym1.fname);
sym1.fname=fullfile(pth,[fn,'_sym',ext]);
ea_write_nii(sym1);

[pth,fn,ext]=fileparts(sym2.fname);
sym2.fname=fullfile(pth,[fn,'_sym',ext]);
ea_write_nii(sym2);

delete(po1)
delete(po2)