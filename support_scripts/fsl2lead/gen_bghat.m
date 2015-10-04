%-----------------------------------------------------------------------
% Job configuration created by cfg_util (rev $Rev: 4252 $)
%-----------------------------------------------------------------------
clear
clc

leaddir=[fileparts(which('lead')),filesep];


        mkdir([leaddir,'atlases',filesep,'BGHAT']);

for region=[1:5]
    switch region
        case 1
            str='rh';
            reg='caudate';
        case 2
            str='rh';
            reg='pallidum';
        case 3
            str='rh';
            reg='GPe';
        case 4
            str='rh';
            reg='GPi';
        case 5
            str='rh';
            reg='STN';
    end
        mkdir([leaddir,'atlases',filesep,'BGHAT',filesep,str]);
        
        matlabbatch{1}.spm.util.imcalc.input = {['BGHATspm.nii,1']};
        matlabbatch{1}.spm.util.imcalc.output = [reg,'.nii'];
        matlabbatch{1}.spm.util.imcalc.outdir = {[leaddir,'atlases',filesep,'BGHAT',filesep,str,filesep]};
        matlabbatch{1}.spm.util.imcalc.expression = ['i1==',num2str(region)];
        matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
        matlabbatch{1}.spm.util.imcalc.options.mask = 0;
        matlabbatch{1}.spm.util.imcalc.options.interp = 3;
        matlabbatch{1}.spm.util.imcalc.options.dtype = 16;
        jobs{1}=matlabbatch;
        cfg_util('run',jobs);
        clear matlabbatch jobs
        
        ea_reslice_nii([leaddir,'atlases',filesep,'BGHAT',filesep,str,filesep,reg,'.nii'],[leaddir,'atlases',filesep,'BGHAT',filesep,str,filesep,reg,'.nii'],[0.5,0.5,0.5],0,0,3)
        
        
                mkdir([leaddir,'atlases',filesep,'BGHAT',filesep,'lh']);

nii=load_nii([leaddir,'atlases',filesep,'BGHAT',filesep,str,filesep,reg,'.nii']);
nii.img=flipdim(nii.img,1);
save_nii(nii,[leaddir,'atlases',filesep,'BGHAT',filesep,'lh',filesep,reg,'.nii']);
   

nii=load_nii([leaddir,'atlases',filesep,'BGHAT',filesep,'lh',filesep,reg,'.nii']);
nii.img=flipdim(nii.img,1);
save_nii(nii,[leaddir,'atlases',filesep,'BGHAT',filesep,'rh',filesep,reg,'.nii']);

end




