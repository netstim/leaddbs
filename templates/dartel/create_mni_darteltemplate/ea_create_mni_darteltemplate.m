function ea_create_mni_darteltemplate(segmentresult)


wd=[fileparts(which('ea_create_mni_darteltemplate')),filesep];
gunzip([wd,'dartelmni_6_hires.nii.gz']);
spm_file_split([wd,'dartelmni_6_hires.nii']);
gs=[0,2,3,5,6,8];
expo=6:-1:1;
for s=1:6

    for tpm=1:3

        % smooth
        if gs(s)

            matlabbatch{1}.spm.spatial.smooth.data = {[wd,'dartelmni_6_hires_',sprintf('%05d',tpm),'.nii,1']};
            matlabbatch{1}.spm.spatial.smooth.fwhm = [gs(s),gs(s),gs(s)];
            matlabbatch{1}.spm.spatial.smooth.dtype = 0;
            matlabbatch{1}.spm.spatial.smooth.im = 0;
            matlabbatch{1}.spm.spatial.smooth.prefix = ['s',num2str(gs(s))];
            jobs{1}=matlabbatch;
            cfg_util('run',jobs);
                    clear jobs matlabbatch

        else

            matlabbatch{1}.spm.util.imcalc.input = {[wd,'dartelmni_6_hires_',sprintf('%05d',tpm),'.nii,1']};
            matlabbatch{1}.spm.util.imcalc.output = [wd,'s0','dartelmni_6_hires_',sprintf('%05d',tpm),'.nii'];
            matlabbatch{1}.spm.util.imcalc.outdir = {wd};
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
        matlabbatch{1}.spm.util.imcalc.input = {[segmentresult,',1'];
            [wd,'s',num2str(gs(s)),'dartelmni_6_hires_',sprintf('%05d',tpm),'.nii']};
        matlabbatch{1}.spm.util.imcalc.output = [wd,'s',num2str(gs(s)),'dartelmni_6_hires_',sprintf('%05d',tpm),'.nii'];
        matlabbatch{1}.spm.util.imcalc.outdir = {wd};
        matlabbatch{1}.spm.util.imcalc.expression = 'i2';
        matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
        matlabbatch{1}.spm.util.imcalc.options.mask = 0;
        matlabbatch{1}.spm.util.imcalc.options.interp = 1;
        matlabbatch{1}.spm.util.imcalc.options.dtype = 16;
        jobs{1}=matlabbatch;
        cfg_util('run',jobs);
        clear jobs matlabbatch



    end

    matlabbatch{1}.spm.util.cat.vols = {[wd,'s',num2str(gs(s)),'dartelmni_6_hires_',sprintf('%05d',1),'.nii'];
        [wd,'s',num2str(gs(s)),'dartelmni_6_hires_',sprintf('%05d',2),'.nii'];
        [wd,'s',num2str(gs(s)),'dartelmni_6_hires_',sprintf('%05d',3),'.nii']
        };
    matlabbatch{1}.spm.util.cat.name = [wd,'dartelmni_',num2str(expo(s)),'.nii'];
    matlabbatch{1}.spm.util.cat.dtype = 0;
    jobs{1}=matlabbatch;
    cfg_util('run',jobs);
    clear jobs matlabbatch

    disp('Cleaning up.');

    % cleanup
    delete([wd,'s',num2str(gs(s)),'dartelmni_6_hires_00*.*']);

end
    % further cleanup
delete([wd,'dartelmni_*.mat']);
delete([wd,'dartelmni_6_hires_',sprintf('%05d',1),'.nii']);
delete([wd,'dartelmni_6_hires_',sprintf('%05d',2),'.nii']);
delete([wd,'dartelmni_6_hires_',sprintf('%05d',3),'.nii']);

gzip([wd,'dartelmni_6_hires.nii']);
delete([wd,'dartelmni_6_hires.nii']);
disp('Done.');

movefile([wd,'dartelmni_1.nii'],[fileparts(fileparts(wd)),filesep,'dartelmni_1.nii'])
movefile([wd,'dartelmni_2.nii'],[fileparts(fileparts(wd)),filesep,'dartelmni_2.nii'])
movefile([wd,'dartelmni_3.nii'],[fileparts(fileparts(wd)),filesep,'dartelmni_3.nii'])
movefile([wd,'dartelmni_4.nii'],[fileparts(fileparts(wd)),filesep,'dartelmni_4.nii'])
movefile([wd,'dartelmni_5.nii'],[fileparts(fileparts(wd)),filesep,'dartelmni_5.nii'])
movefile([wd,'dartelmni_6.nii'],[fileparts(fileparts(wd)),filesep,'dartelmni_6.nii'])


