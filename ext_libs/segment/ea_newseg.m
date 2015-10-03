function ea_newseg(directory,file,dartel,options)


disp('Segmentation...');
    switch spm('ver')
        
        case 'SPM8'
           
            tpminf=fullfile(fileparts(which('spm')),'toolbox','Seg',['TPM.nii']);
            tpmoutf=[options.earoot,'templates',filesep,'TPM.nii'];
            if ~exist(tpmoutf,'file')
                ea_reslice_nii(tpminf,tpmoutf,[0.5,0.5,0.5],3);
            end
            
            matlabbatch{1}.spm.tools.preproc8.channel.vols = {[directory,file,',1']};
            matlabbatch{1}.spm.tools.preproc8.channel.biasreg = 0.0001;
            matlabbatch{1}.spm.tools.preproc8.channel.biasfwhm = 60;
            matlabbatch{1}.spm.tools.preproc8.channel.write = [0 0];
            for tpm=1:6
                
                matlabbatch{1}.spm.tools.preproc8.tissue(tpm).tpm = {[tpmoutf,',',num2str(tpm)]};
                matlabbatch{1}.spm.tools.preproc8.tissue(tpm).ngaus = 2;
                matlabbatch{1}.spm.tools.preproc8.tissue(tpm).native = [1 0];
                matlabbatch{1}.spm.tools.preproc8.tissue(tpm).warped = [0 0];
                
                if dartel && tpm < 4
                matlabbatch{1}.spm.tools.preproc8.tissue(tpm).native = [1 1];
                end
                
            end
            
            matlabbatch{1}.spm.tools.preproc8.warp.mrf = 0;
            matlabbatch{1}.spm.tools.preproc8.warp.reg = 4;
            matlabbatch{1}.spm.tools.preproc8.warp.affreg = 'mni';
            matlabbatch{1}.spm.tools.preproc8.warp.samp = 3;
            matlabbatch{1}.spm.tools.preproc8.warp.write = [1 1];
            cfg_util('run',{matlabbatch});
            clear matlabbatch
        case 'SPM12'
                load([options.earoot,'ext_libs',filesep,'segment',filesep,'segjob12']);
            tpminf=fullfile(fileparts(which('spm')),'tpm',['TPM.nii']);
            job.channel.vols{1}=[directory,file,',1'];
            for tpm=1:6
                job.tissue(tpm).tpm=[tpminf,',',num2str(tpm)]; 
                if dartel && tpm < 4
                job.tissue(tpm).native(2)=1;
                end
            end
                        job.warp.write=[1,1];
            spm_preproc_run(job); % run "Segment" in SPM 12 (Old "Segment" is now referred to as "Old Segment").
    end

    disp('Done.');