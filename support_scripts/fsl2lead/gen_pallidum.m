clear
clc

% put the following files (unzipped versions) from your FSL installation
% inside the folder of this script and run it:
% HarvardOxford-cort-maxprob-thr0-1mm.nii
% HarvardOxford-sub-maxprob-thr0-1mm.nii


leaddir=[fileparts(which('lead')),filesep];
mkdir([leaddir,'atlases',filesep,'Pallidum']);

for lr=[13,52] % pallidum
    switch lr
        case 13
            str='lh';
        case 52
            str='rh';
    end
    for t=[0]
        
        matlabbatch{1}.spm.util.imcalc.input = {['HarvardOxford-sub-maxprob-thr',num2str(t),'-1mm.nii,1']};
        matlabbatch{1}.spm.util.imcalc.output = ['pall_',num2str(t),'.nii'];
        matlabbatch{1}.spm.util.imcalc.outdir = {[leaddir,'atlases',filesep,'Pallidum',filesep,str,filesep]};
        matlabbatch{1}.spm.util.imcalc.expression = ['i1==',num2str(lr)];
        matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
        matlabbatch{1}.spm.util.imcalc.options.mask = 0;
        matlabbatch{1}.spm.util.imcalc.options.interp = 0;
        matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
        jobs{1}=matlabbatch;
        cfg_util('run',jobs);
        clear matlabbatch jobs        
        
        
    end
end