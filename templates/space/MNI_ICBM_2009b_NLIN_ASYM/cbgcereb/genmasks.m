%ea_conformspaceto('TPM_00001.nii','atropos_thal.nii',1);
%ea_conformspaceto('TPM_00001.nii','atropos_cb.nii',1);
%ea_conformspaceto('TPM_00001.nii','atropos_ct.nii',1);
%ea_conformspaceto('TPM_00001.nii','atropos_bs.nii',1);

c1=ea_load_nii('TPM_00001.nii');
cb=ea_load_nii('atropos_cb.nii');
ct=ea_load_nii('atropos_ct.nii');
bs=ea_load_nii('atropos_bs.nii');
thal=ea_load_nii('atropos_thal.nii');


c1cb=c1;
c1cb.img=c1cb.img.*double(cb.img);
c1cb.fname='cerebellum.nii';
ea_write_nii(c1cb);

c1cb=c1;
c1cb.img=c1cb.img.*double(bs.img);
c1cb.fname='brainstem.nii';
ea_write_nii(c1cb);

c1cb=c1;
c1cb.img=c1cb.img.*double(ct.img);
c1cb.fname='cortex.nii';
ea_write_nii(c1cb);

c1cb=c1;
c1cb.img=c1cb.img.*double(thal.img);
c1cb.fname='basal_ganlia.nii';
ea_write_nii(c1cb);
