darteltemp=''; % define path to dartel template matching your data.

gs=[0,2,3,5,6,8];
for s=1:6
    
    for tpm=1:3
        
        % smooth
        if gs(s)
            matlabbatch{1}.spm.spatial.smooth.data = {['dartelmni_6_hires.nii,',num2str(tpm)]};
            matlabbatch{1}.spm.spatial.smooth.fwhm = [gs(s),gs(s),gs(s)];
            matlabbatch{1}.spm.spatial.smooth.dtype = 0;
            matlabbatch{1}.spm.spatial.smooth.im = 0;
            matlabbatch{1}.spm.spatial.smooth.prefix = ['s',num2str(gs(s)),'c',num2str(tpm)];
            jobs{1}=matlabbatch;
            cfg_util('run',jobs);
        else
            matlabbatch{1}.spm.util.imcalc.input = {['dartelmni_6_hires.nii,',num2str(tpm)]};
            matlabbatch{1}.spm.util.imcalc.output = ['s0','c',num2str(tpm),'dartelmni_6_hires.nii'];
            matlabbatch{1}.spm.util.imcalc.outdir = {''};
            matlabbatch{1}.spm.util.imcalc.expression = 'i1';
            matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
            matlabbatch{1}.spm.util.imcalc.options.mask = 0;
            matlabbatch{1}.spm.util.imcalc.options.interp = 1;
            matlabbatch{1}.spm.util.imcalc.options.dtype = 32;
            jobs{1}=matlabbatch;
            cfg_util('run',jobs);
        end
        clear jobs matlabbach
        
        % set to resolution of darteltemp file
        matlabbatch{1}.spm.util.imcalc.input = {darteltemp
            ['s',num2str(gs(s)),'c',num2str(tpm),'dartelmni_6_hires.nii']};
        matlabbatch{1}.spm.util.imcalc.output = ['s',num2str(gs(s)),'c',num2str(tpm),'dartelmni_6_hires.nii'];
        matlabbatch{1}.spm.util.imcalc.outdir = {''};
        matlabbatch{1}.spm.util.imcalc.expression = 'i2';
        matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
        matlabbatch{1}.spm.util.imcalc.options.mask = 0;
        matlabbatch{1}.spm.util.imcalc.options.interp = 1;
        matlabbatch{1}.spm.util.imcalc.options.dtype = 32;
        jobs{1}=matlabbatch;
        cfg_util('run',jobs);
        clear jobs matlabbach
        
        
    end
    
    matlabbatch{1}.spm.util.cat.vols = {['s',num2str(gs(s)),'c',num2str(1),'dartelmni_6_hires.nii']
        ['s',num2str(gs(s)),'c',num2str(2),'dartelmni_6_hires.nii']
        ['s',num2str(gs(s)),'c',num2str(3),'dartelmni_6_hires.nii']
        };
    matlabbatch{1}.spm.util.cat.name = ['dartelmni_',num2str(tpm),'.nii'];
    matlabbatch{1}.spm.util.cat.dtype = 0;
    jobs{1}=matlabbatch;
    cfg_util('run',jobs);
    clear jobs matlabbach
    
end