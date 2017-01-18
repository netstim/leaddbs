function ea_redo_inv(varargin)

directory=varargin{1};
options=varargin{2};

if nargin>2
    forinv=varargin{3};
else
    forinv='inverse';
end

switch forinv
    case 'inverse'
        matlabbatch{1}.spm.util.defs.comp{1}.inv.comp{1}.def = {[directory,'y_ea_normparams.nii']};
        matlabbatch{1}.spm.util.defs.comp{1}.inv.space = {[directory,options.prefs.prenii_unnormalized]};
        matlabbatch{1}.spm.util.defs.out{1}.savedef.ofname = 'ea_inv_normparams.nii';
        matlabbatch{1}.spm.util.defs.out{1}.savedef.savedir.saveusr = {directory};
        spm_jobman('run',{matlabbatch});
        
    case 'forward'    
        matlabbatch{1}.spm.util.defs.comp{1}.id.space = {[ea_space(options),options.primarytemplate,'.nii']};
        matlabbatch{1}.spm.util.defs.comp{1}.def = {[directory,'y_ea_normparams.nii']};
        matlabbatch{1}.spm.util.defs.out{1}.savedef.ofname = 'ea_normparams.nii';
        matlabbatch{1}.spm.util.defs.out{1}.savedef.savedir.saveusr = {directory};
        spm_jobman('run',{matlabbatch});
     
end
