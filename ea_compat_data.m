function ea_compat_data

if exist([ea_getearoot,'templates',filesep,'mni_hires_t1.nii'],'file');
    movefile([ea_getearoot,'templates'], [ea_getearoot,'templates_temp']);
    mkdir([ea_getearoot,'templates',filesep,'space',filesep,'MNI_ICBM_2009b_NLIN_ASYM']);
    movefile([ea_getearoot,'templates_temp'],[ea_getearoot,'templates',filesep,'space',filesep,'MNI_ICBM_2009b_NLIN_ASYM',filesep,'templates']);
    
    movefile([ea_getearoot,'atlases'],[ea_getearoot,'templates',filesep,'space',filesep,'MNI_ICBM_2009b_NLIN_ASYM',filesep,'atlases']);
    
    movefile([ea_space,'mni_hires_t1.nii'],[ea_space,'t1.nii']);
    movefile([ea_space,'mni_hires_t2.nii'],[ea_space,'t2.nii']);
    movefile([ea_space,'mni_hires_pd.nii'],[ea_space,'pd.nii']);
    movefile([ea_space,'mni_hires_fa.nii'],[ea_space,'fa.nii']);
    movefile([ea_space,'mni_hires_bb.nii'],[ea_space,'bb.nii']);
    movefile([ea_space,'TPM_2009b.nii'],[ea_space,'TPM.nii']);
    movefile([ea_space,'mni_hires_distal.nii'],[ea_space,'distal.nii']);
    movefile([ea_space,'mni_hires_wires.mat'],[ea_space,'wires.mat']);
end


