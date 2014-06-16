%-----------------------------------------------------------------------
% Job configuration created by cfg_util (rev $Rev: 4252 $)
%-----------------------------------------------------------------------
root='/PA/Neuro/_projects/lead/lead/templates/segment/';
ages={'E','M','Y'};
sides={'R','L'};
for age=1:length(ages)
    for side=1:length(sides);
        nii=load_untouch_nii([root,'STN_',sides{side},'_',ages{age},'.nii']);
        nii.img=double(nii.img);
        nii.img=nii.img/max(nii.img(:)); % max 1
        save_untouch_nii(nii,[root,'eSTN_',sides{side},'_',ages{age},'.nii']);
        clear nii
        
        matlabbatch{1}.spm.util.imcalc.input = {
            [root,'TPM.nii,1']
            [root,'eSTN_',sides{side},'_',ages{age},'.nii,1']
            };
        matlabbatch{1}.spm.util.imcalc.output = ['tSTN_',sides{side},'_',ages{age},'.nii'];
        matlabbatch{1}.spm.util.imcalc.outdir = {root};
        matlabbatch{1}.spm.util.imcalc.expression = 'i2';
        matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
        matlabbatch{1}.spm.util.imcalc.options.mask = 0;
        matlabbatch{1}.spm.util.imcalc.options.interp = -4;
        matlabbatch{1}.spm.util.imcalc.options.dtype = 16;
        
        
        jobs{1}=matlabbatch;
        cfg_util('run',jobs);
        clear jobs matlabbatch
        
        
        
    end
    tnii=load_untouch_nii([root,'TPM.nii']);
    rnii=load_untouch_nii([root,'tSTN_',sides{1},'_',ages{age},'.nii']);
    lnii=load_untouch_nii([root,'tSTN_',sides{2},'_',ages{age},'.nii']);
    tnii.img(:,:,:,1)=tnii.img(:,:,:,1)+lnii.img+rnii.img;
    save_untouch_nii(tnii,[root,'STN_',ages{age},'_TPM.nii']);
    
    
end