function ea_unpackspace

disp(['Unpacking space ',ea_getspace,'...']);
disp('This could take a while...');
try gunzip([ea_space,'*.nii.gz']); end
try delete([ea_space,'*.nii.gz']); end
load([ea_space,'spacedef.mat']);
ea_genauxspace([],[],spacedef.defaultatlas);

ea_delete([ea_space,'need_build']);

disp('Unpacking done.');
