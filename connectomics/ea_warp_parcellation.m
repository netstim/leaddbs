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

            ea_ants_applytransforms(options, ...
                {[ea_space(options,'labeling'),options.lc.general.parcellation,'.nii']}, ...
                {[options.root,options.patientname,filesep,'templates',filesep,'labeling',filesep,'w',options.lc.general.parcellation,'.nii']},...
                1,'','','GenericLabel');

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
    copyfile([options.root,options.patientname,filesep,options.prefs.prenii_unnormalized],[options.root,options.patientname,filesep,'c',options.prefs.prenii_unnormalized]);
    matlabbatch{1}.spm.spatial.coreg.estwrite.ref = {[options.root,options.patientname,filesep,ref_filename,',1']};
    matlabbatch{1}.spm.spatial.coreg.estwrite.source = {[options.root,options.patientname,filesep,'c',options.prefs.prenii_unnormalized,',1']};
    matlabbatch{1}.spm.spatial.coreg.estwrite.other = {[directory,'templates',filesep,'labeling',filesep,'w',options.lc.general.parcellation,'.nii,1']};
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = 0;
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.mask = 0;
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.prefix = ['r',b0rest];
    spm_jobman('run',{matlabbatch});
    clear matlabbatch

    delete([options.root,options.patientname,filesep,'c',options.prefs.prenii_unnormalized]);
    delete([options.root,options.patientname,filesep,'r',b0rest,'c',options.prefs.prenii_unnormalized]);
end
