function ea_unpackspace

disp(['Unpacking space ',ea_getspace,'...']);
disp('This could take a while...');
gunzip([ea_space,'*.nii.gz']);
delete([ea_space,'*.nii.gz']);
load([ea_space,'ea_space_def.mat']);
ea_genauxspace([],[],templates.defaultatlas);

delete([ea_space,'packed']);

disp('Unpacking done.');

