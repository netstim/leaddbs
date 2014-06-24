clear
clc
darteltemp='Template_6.nii'; % define path to dartel template matching your data.
spm_file_split('dartelmni_6_hires.nii');
gs=[0,2,3,5,6,8];
expo=6:-1:1;
for s=1:6
    
    for tpm=1:3
        
        % smooth
        if gs(s)
            
            matlabbatch{1}.spm.spatial.smooth.data = {['dartelmni_6_hires_',sprintf('%05d',tpm),'.nii,1']};
            matlabbatch{1}.spm.spatial.smooth.fwhm = [gs(s),gs(s),gs(s)];
            matlabbatch{1}.spm.spatial.smooth.dtype = 0;
            matlabbatch{1}.spm.spatial.smooth.im = 0;
            matlabbatch{1}.spm.spatial.smooth.prefix = ['s',num2str(gs(s))];
            jobs{1}=matlabbatch;
            cfg_util('run',jobs);
                    clear jobs matlabbatch

        else
            
            matlabbatch{1}.spm.util.imcalc.input = {['dartelmni_6_hires_',sprintf('%05d',tpm),'.nii,1']};
            matlabbatch{1}.spm.util.imcalc.output = ['s0','dartelmni_6_hires_',sprintf('%05d',tpm),'.nii'];
            matlabbatch{1}.spm.util.imcalc.outdir = {''};
            matlabbatch{1}.spm.util.imcalc.expression = 'i1';
            matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
            matlabbatch{1}.spm.util.imcalc.options.mask = 0;
            matlabbatch{1}.spm.util.imcalc.options.interp = 1;
            matlabbatch{1}.spm.util.imcalc.options.dtype = 16;
            jobs{1}=matlabbatch;
            cfg_util('run',jobs);
                    clear jobs matlabbatch

        end
        clear jobs matlabbatch
        
        % set to resolution of darteltemp file
        matlabbatch{1}.spm.util.imcalc.input = {[darteltemp,',1']
            ['s',num2str(gs(s)),'dartelmni_6_hires_',sprintf('%05d',tpm),'.nii']};
        matlabbatch{1}.spm.util.imcalc.output = ['s',num2str(gs(s)),'dartelmni_6_hires_',sprintf('%05d',tpm),'.nii'];
        matlabbatch{1}.spm.util.imcalc.outdir = {''};
        matlabbatch{1}.spm.util.imcalc.expression = 'i2';
        matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
        matlabbatch{1}.spm.util.imcalc.options.mask = 0;
        matlabbatch{1}.spm.util.imcalc.options.interp = 1;
        matlabbatch{1}.spm.util.imcalc.options.dtype = 16;
        jobs{1}=matlabbatch;
        cfg_util('run',jobs);
        clear jobs matlabbatch
        
        
        
    end
    
    matlabbatch{1}.spm.util.cat.vols = {['s',num2str(gs(s)),'dartelmni_6_hires_',sprintf('%05d',1),'.nii']
        ['s',num2str(gs(s)),'dartelmni_6_hires_',sprintf('%05d',2),'.nii']
        ['s',num2str(gs(s)),'dartelmni_6_hires_',sprintf('%05d',3),'.nii']
        };
    matlabbatch{1}.spm.util.cat.name = ['dartelmni_',num2str(expo(s)),'.nii'];
    matlabbatch{1}.spm.util.cat.dtype = 0;
    jobs{1}=matlabbatch;
    cfg_util('run',jobs);
    clear jobs matlabbatch
    % cleanup
    delete(['s',num2str(gs(s)),'dartelmni_6_hires_00*.*']);

end
    % further cleanup
delete('dartelmni_*.mat');
delete(['dartelmni_6_hires_',sprintf('%05d',1),'.nii']);
delete(['dartelmni_6_hires_',sprintf('%05d',2),'.nii']);   
delete(['dartelmni_6_hires_',sprintf('%05d',3),'.nii']);