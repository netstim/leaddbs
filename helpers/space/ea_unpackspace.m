function ea_unpackspace

disp(['Unpacking space ',ea_getspace,'...']);
disp('This could take a while...');
gunzip([ea_space,'*.nii.gz']);
delete([ea_space,'*.nii.gz']);
ea_genauxspace([],[],ea_getspace);

delete([ea_space,'packed']);

disp('Unpacking done.');

