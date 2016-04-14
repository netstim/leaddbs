function doseed2map_lite(root)
% sample batch that will generate a connectivity map using any *_seed.nii
% file as seed and the HCP_MGH_30fold group connectome as data basis.



if ~isdeployed
addpath(genpath('/autofs/cluster/nimlab/USERS/Andy/lead_dbs'));
addpath(genpath('/autofs/cluster/nimlab/USERS/Andy/spm12'));
addpath(genpath('/autofs/cluster/nimlab/USERS/Andy/MNI_connectomes'));
end
if ~strcmp(root(end),filesep)
    root=[root,filesep];
end
fdir=dir([root,'*_seed.nii']);

for fi=1:length(fdir)
   fis{fi}=[root,fdir(fi).name]; 
end

ea_seed2map_weighted('/autofs/cluster/nimlab/USERS/Andy/MNI_connectomes/HCP_MGH_30fold_groupconnectome_gqi_lite.mat',...
    fis,...
    '/autofs/cluster/nimlab/USERS/Andy/MNI_connectomes/template_to_use.nii');


exit
