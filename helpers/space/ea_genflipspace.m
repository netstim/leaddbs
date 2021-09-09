function ea_genflipspace
if isempty(ea_regexpdir([ea_space,'fliplr'],'glanatComposite\.(h5|nii\.gz)',0))
    
    answ=questdlg('The transform to map from left to right hemisphere in 2009b NLIN ASYM space needs to be installed or built. For reproducable results please install the transform',...
    'Install Flip Transform','Install','Build','Cancel','Install');
    switch answ
        case 'Cancel'
            ea_error('User pressed cancel.');
        case 'Install'
            ea_checkinstall('nlinflip');
            return
        otherwise
            ea_dispt('Generating nonlinear warp between hemispheres. This may take a while...');
            mkdir([ea_space,'fliplr']);
            spacedef=ea_getspacedef;
            for t=1:length(spacedef.templates)
                copyfile([ea_space,spacedef.templates{t},'.nii'],[ea_space,'fliplr',filesep,'anat_',spacedef.templates{t},'.nii']);
                ea_flip_lr([ea_space,'fliplr',filesep,'anat_',spacedef.templates{t},'.nii'],[ea_space,'fliplr',filesep,'anat_',spacedef.templates{t},'.nii']);
            end
            
            options=ea_getptopts([ea_space,'fliplr',filesep]);
            options.coregmr.method='Do not coregister MRIs (already coregistered)';
            options.modality=1;
            options.overwriteapproved=1;
            ea_dumpmethod(options, 'norm');
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

end
