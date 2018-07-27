function ea_warp_parcellation(ref_filename,b0rest,options)
directory=[options.root,options.patientname,filesep];

if ~exist([directory,'templates',filesep,'labeling',filesep,'w',options.lc.general.parcellation,'.nii'],'file')
    %% warp atlas into pre_tra-space:
    if ~exist([directory,'templates'],'dir')
        mkdir([directory,'templates']);
    end
    if ~exist([directory,'templates',filesep,'labeling'],'dir')
        mkdir([directory,'templates',filesep,'labeling']);
    end

    whichnormmethod=ea_whichnormmethod([options.root,options.patientname,filesep]);
    switch whichnormmethod
        case ea_getantsnormfuns
            useinterp='GenericLabel';
            parc=ea_load_nii([ea_space(options,'labeling'),options.lc.general.parcellation,'.nii']);
            if length(unique(parc.img(:)))>500
                useinterp='NearestNeighbor'; % both GenericLabel and MultiLabel take ages on high dimensional parcellations.
            end
            ea_ants_applytransforms(options, ...
                {[ea_space(options,'labeling'),options.lc.general.parcellation,'.nii']}, ...
                {[options.root,options.patientname,filesep,'templates',filesep,'labeling',filesep,'w',options.lc.general.parcellation,'.nii']},...
                1,'','',useinterp);

        case ea_getfslnormfuns

            ea_fsl_applytransforms(options, ...
                {[ea_space(options,'labeling'),options.lc.general.parcellation,'.nii']}, ...
                {[options.root,options.patientname,filesep,'templates',filesep,'labeling',filesep,'w',options.lc.general.parcellation,'.nii']},...
                1,'','','nn');

        otherwise
            switch spm('ver')
                case 'SPM8'


                    matlabbatch{1}.spm.util.defs.comp{1}.def = {[options.root,options.patientname,filesep,'y_ea_inv_normparams.nii']};
                    matlabbatch{1}.spm.util.defs.ofname = '';
                    matlabbatch{1}.spm.util.defs.fnames = {[ea_space(options,'labeling'),options.lc.general.parcellation,'.nii,1']};
                    matlabbatch{1}.spm.util.defs.savedir.saveusr = {[options.root,options.patientname,filesep,'templates',filesep,'labeling',filesep]};
                    matlabbatch{1}.spm.util.defs.interp = 0;
                    spm_jobman('run',{matlabbatch});
                    clear matlabbatch

                case 'SPM12'

                    matlabbatch{1}.spm.util.defs.comp{1}.def = {[options.root,options.patientname,filesep,'y_ea_inv_normparams.nii']};
                    matlabbatch{1}.spm.util.defs.out{1}.pull.fnames = {[ea_space(options,'labeling'),options.lc.general.parcellation,'.nii']};
                    matlabbatch{1}.spm.util.defs.out{1}.pull.savedir.saveusr = {[options.root,options.patientname,filesep,'templates',filesep,'labeling',filesep]};
                    matlabbatch{1}.spm.util.defs.out{1}.pull.interp = 0;
                    matlabbatch{1}.spm.util.defs.out{1}.pull.mask = 1;
                    matlabbatch{1}.spm.util.defs.out{1}.pull.fwhm = [0 0 0];
                    spm_jobman('run',{matlabbatch});
                    clear matlabbatch

            end
    end
end

if ~exist([directory,'templates',filesep,'labeling',filesep,'r',b0rest,'w',options.lc.general.parcellation,'.nii'],'file')
    %% coreg atlas into b0-space:
    [~,reffn]=fileparts(ref_filename);
    copyfile([options.root,options.patientname,filesep,options.prefs.prenii_unnormalized],[options.root,options.patientname,filesep,'c',options.prefs.prenii_unnormalized]);
    copyfile([directory,'templates',filesep,'labeling',filesep,'w',options.lc.general.parcellation,'.nii'],...
        [directory,'templates',filesep,'labeling',filesep,'r',reffn,'w',options.lc.general.parcellation,'.nii']);

    ea_coreg2images(options,[options.root,options.patientname,filesep,'c',options.prefs.prenii_unnormalized],...
        [options.root,options.patientname,filesep,ref_filename],...
        [options.root,options.patientname,filesep,'c',options.prefs.prenii_unnormalized],...
        {[directory,'templates',filesep,'labeling',filesep,'r',reffn,'w',options.lc.general.parcellation,'.nii']},0,[],0);


    ea_gencheckregpair([directory,'templates',filesep,'labeling',filesep,'r',reffn,'w',options.lc.general.parcellation],...
        [options.root,options.patientname,filesep,reffn],...
        [options.root,options.patientname,filesep,'checkreg',filesep,options.lc.general.parcellation,'2',reffn,'.png']);

    ea_delete([options.root,options.patientname,filesep,'c',options.prefs.prenii_unnormalized]);
    ea_delete([options.root,options.patientname,filesep,'r',reffn,'c',options.prefs.prenii_unnormalized]);
end
