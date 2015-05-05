function ea_newseg(directory,file,dartel,options)


disp('Segmentation...');
    switch spm('ver')
        
        case 'SPM8'
            load([options.earoot,'ext_libs',filesep,'segment',filesep,'segjob8']);
            % make some reslicing to obtain higher resolution in SPM 8.
            job.channel.vols{1}=[directory,file,',1'];
            tpminf=fullfile(fileparts(which('spm')),'toolbox','Seg',['TPM.nii']);
            tpmoutf=[options.earoot,'templates',filesep,'TPM.nii'];
            if ~exist(tpmoutf,'file')
                reslice_nii(tpminf,tpmoutf,[0.5,0.5,0.5],3);
            end
            for tpm=1:6
                job.tissue(tpm).tpm=[tpmoutf,',',num2str(tpm)];
                if dartel && tpm < 4
                job.tissue(tpm).native(2)=1;
                end                    
                
            end
            job.resolution=0.5;
            job.warp.write=[1,1];
            ea_spm_preproc_run(job); % exactly the same as the SPM version ("New Segment" in SPM8) but with an increase in resolution to 0.5 mm iso.
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