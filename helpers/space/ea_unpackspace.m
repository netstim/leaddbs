function ea_unpackspace

gunzip([ea_space,'*.nii.gz']);
delete([ea_space,'*.nii.gz']);
            ea_genauxspace([],[],ea_getspace);

            delete([ea_space,'packed']);