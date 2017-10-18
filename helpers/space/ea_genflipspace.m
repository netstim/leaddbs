function ea_genflipspace
if ~exist([ea_space,'fliplr',filesep,'glanatComposite.h5'],'file');
    ea_dispt('Generating nonlinear warp between hemispheres');
    mkdir([ea_space,'fliplr']);
    spacedef=ea_getspacedef;
    for t=1:length(spacedef.templates)
        copyfile([ea_space,spacedef.templates{t},'.nii'],[ea_space,'fliplr',filesep,'anat_',spacedef.templates{t},'.nii']);
        ea_flip_lr([ea_space,'fliplr',filesep,'anat_',spacedef.templates{t},'.nii'],[ea_space,'fliplr',filesep,'anat_',spacedef.templates{t},'.nii']);
        
    end
    
    options=ea_getptopts([ea_space,'fliplr',filesep]);
    options.coregmr.method='Do not coregister MRIs (already coregistered)';
    options.modality=1;
    ea_dumpnormmethod(options,'ea_normalize_ants','normmethod'); % has to come first due to applynormalization.
    ea_normalize_ants(options,0);
    
    for t=1:length(spacedef.templates)
        delete([ea_space,'fliplr',filesep,'anat_',spacedef.templates{t},'.nii']);
    end
    delete([ea_space,'fliplr',filesep,'glanat.nii']);
    delete([ea_space,'fliplr',filesep,'lanat.nii']);
    delete([ea_space,'fliplr',filesep,'ea_methods.txt']);
    delete([ea_space,'fliplr',filesep,'ea_ants_command.txt']);
    ea_dispt('');
end