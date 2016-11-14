function   ea_apply_normalization_tofile(options,from,to,directory,useinverse)
% this function applies lead-dbs normalizations to nifti files.

switch ea_whichnormmethod(directory);
    case ea_getantsnormfuns % ANTs part here
        
        ea_ants_applytransforms(options,from,to,useinverse);
        
    case ea_getfslnormfuns % FSL part here
        
        ea_fsl_applytransforms(options,from,to,useinverse);
        
    otherwise % SPM part here
        for fi=1:length(from) % assume from and to have same length (must have for this to work)
            matlabbatch{1}.spm.util.defs.comp{1}.def = {[directory,'y_ea_normparams.nii']}; % non-inverse usage okay here since using pushforward method.
            matlabbatch{1}.spm.util.defs.out{1}.push.fnames = from(fi);
            matlabbatch{1}.spm.util.defs.out{1}.push.weight = {''};
            matlabbatch{1}.spm.util.defs.out{1}.push.savedir.saveusr = {fileparts(to{fi})};
            matlabbatch{1}.spm.util.defs.out{1}.push.fov.file = {[directory,options.prefs.prenii_unnormalized]};
            matlabbatch{1}.spm.util.defs.out{1}.push.preserve = 0;
            matlabbatch{1}.spm.util.defs.out{1}.push.fwhm = [0 0 0];
            matlabbatch{1}.spm.util.defs.out{1}.push.prefix = '';
            spm_jobman('run',{matlabbatch});
            clear matlabbatch
            [pth]=fileparts(to{fi});
            [~,fn,ext]=fileparts(from{fi});
            movefile(fullfile(pth,['w',fn,ext]),to{fi});
        end
        
end