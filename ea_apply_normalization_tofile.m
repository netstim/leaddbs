function ea_apply_normalization_tofile(options,from,to,directory,useinverse,interp,refim)
% this function applies lead-dbs normalizations to nifti files.
% currently just used to generate patient specific atlases,i.e., from MNI
% space to native space

if ~exist('interp','var')
    interp=4;
end

if ~exist('refim','var')
    refim='';
end
switch ea_whichnormmethod(directory)
    case ea_getantsnormfuns % ANTs part here
        
        ea_ants_applytransforms(options,from,to,useinverse,refim,'',interp);
        
    case ea_getfslnormfuns % FSL part here
        
        ea_fsl_applytransforms(options,from,to,useinverse);
        
    otherwise % SPM part here
        for fi=1:length(from) % assume from and to have same length (must have for this to work)
            if strcmp(from{fi}(end-2:end),'.gz') % .gz support
                [pth,fn,ext]=fileparts(from{fi});
                copyfile(from{fi},[tempdir,fn,'.gz']);
                gunzip([tempdir,fn,'.gz']);
                from{fi}=[tempdir,fn];
                wasgz=1;
            else
                wasgz=0;
            end
            if useinverse
               
                matlabbatch{1}.spm.util.defs.comp{1}.def = {[directory,'y_ea_inv_normparams.nii']};
                
                matlabbatch{1}.spm.util.defs.out{1}.pull.fnames = from(fi);
                matlabbatch{1}.spm.util.defs.out{1}.pull.savedir.saveusr = {fileparts(to{fi})};
                matlabbatch{1}.spm.util.defs.out{1}.pull.interp = interp;
                matlabbatch{1}.spm.util.defs.out{1}.pull.mask = 1;
                matlabbatch{1}.spm.util.defs.out{1}.pull.fwhm = [0 0 0];
                matlabbatch{1}.spm.util.defs.out{1}.pull.prefix = '';
                
                
                
                
                             
                
                
                
                
                
                
                
                spm_jobman('run',{matlabbatch});
                clear matlabbatch
                [pth]=fileparts(to{fi});
                [~,fn,ext]=fileparts(from{fi});
                try % fails if to is a w prefixed file already
                    movefile(fullfile(pth,['w',fn,ext]),to{fi});
                end
            else
                
                matlabbatch{1}.spm.util.defs.comp{1}.def = {[directory,'y_ea_normparams.nii']}; 
                matlabbatch{1}.spm.util.defs.out{1}.pull.fnames = from(fi);
                matlabbatch{1}.spm.util.defs.out{1}.pull.savedir.saveusr = {fileparts(to{fi})};
                matlabbatch{1}.spm.util.defs.out{1}.pull.interp = interp;
                matlabbatch{1}.spm.util.defs.out{1}.pull.mask = 1;
                matlabbatch{1}.spm.util.defs.out{1}.pull.fwhm = [0 0 0];
                matlabbatch{1}.spm.util.defs.out{1}.pull.prefix = '';
                
 
                spm_jobman('run',{matlabbatch});
                clear matlabbatch
                [pth]=fileparts(to{fi});
                [~,fn,ext]=fileparts(from{fi});
                try % fails if to is a w prefixed file already
                    movefile(fullfile(pth,['w',fn,ext]),to{fi});
                end
            end
            if wasgz
               gzip(to{fi});
               delete(to{fi});
            end
        end       
end
