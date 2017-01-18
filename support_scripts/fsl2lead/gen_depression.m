%-----------------------------------------------------------------------
% Job configuration created by cfg_util (rev $Rev: 4252 $)
%-----------------------------------------------------------------------
clear
clc


% put the following files (unzipped versions) from your FSL installation
% inside the folder of this script and run it:
% HarvardOxford-cort-maxprob-thr0-1mm.nii
% HarvardOxford-sub-maxprob-thr0-1mm.nii

mkdir([ea_space([],'atlases'),'Depression']);

for region=[27,29,28,33]
    switch region
        case 27
            str='midline';
            reg='SCC';
        case 29
            str='midline';
            reg='aACC';
            case 28
            str='midline';
            reg='PaC';
        case 33
            str='midline';

            reg='FOC';
    end
        mkdir([ea_space([],'atlases'),'Depression',filesep,str]);

        matlabbatch{1}.spm.util.imcalc.input = {['HarvardOxford-cort-maxprob-thr50-1mm.nii,1']};
        matlabbatch{1}.spm.util.imcalc.output = [reg,'.nii'];
        matlabbatch{1}.spm.util.imcalc.outdir = {[ea_space([],'atlases'),'Depression',filesep,str,filesep]};
        matlabbatch{1}.spm.util.imcalc.expression = ['i1==',num2str(region)];
        matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
        matlabbatch{1}.spm.util.imcalc.options.mask = 0;
        matlabbatch{1}.spm.util.imcalc.options.interp = 0;
        matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
        jobs{1}=matlabbatch;
        spm_jobman('run',jobs);
        clear matlabbatch jobs



end


%% add n. acc


for lr=[26,58] % accumbens
    switch lr
        case 26
            str='lh';
        case 58
            str='rh';
    end

        mkdir([ea_space([],'atlases'),'Depression',filesep,str]);
        matlabbatch{1}.spm.util.imcalc.input = {['HarvardOxford-sub-maxprob-thr50-1mm.nii,1']};
        matlabbatch{1}.spm.util.imcalc.output = ['N_Acc','.nii'];
        matlabbatch{1}.spm.util.imcalc.outdir = {[ea_space([],'atlases'),'Depression',filesep,str,filesep]};
        matlabbatch{1}.spm.util.imcalc.expression = ['i1==',num2str(lr)];
        matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
        matlabbatch{1}.spm.util.imcalc.options.mask = 0;
        matlabbatch{1}.spm.util.imcalc.options.interp = 0;
        matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
        jobs{1}=matlabbatch;
        spm_jobman('run',jobs);
        clear matlabbatch jobs





end


