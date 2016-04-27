function doseed2map_weighted(varargin)
% sample batch that will generate a connectivity map using any *_seed.nii
% file as seed and the HCP_MGH_30fold group connectome as data basis.

root=varargin{1};

if nargin>1
    maxdist=varargin{2};
else
    maxdist=[]; % use voxel size later.
end

if ischar(maxdist)
    maxdist=str2double(maxdist);
end

if ~isdeployed
addpath(genpath('/autofs/cluster/nimlab/USERS/Andy/lead_dbs'));
addpath(genpath('/autofs/cluster/nimlab/USERS/Andy/spm12'));
addpath(genpath('/autofs/cluster/nimlab/USERS/Andy/MNI_connectomes'));
end

if strcmp(root(end-2:end),'.gz')
    gunzip(root);
    root=root(1:end-3);
    fis{1}=root;
elseif strcmp(root(end-3:end),'.nii')
    fis{1}=root;
elseif strcmp(root(end-3:end),'.txt')
    fID=fopen(root);
    fis=textscan(fID,'%s');
else % folder supplied, search for _seed.nii files.
    
    if ~strcmp(root(end),filesep)
        root=[root,filesep];
    end
    
    try gunzip([root,'*.nii.gz']); end
    
    fdir=dir([root,'*.nii']);
    
    for fi=1:length(fdir)
        fis{fi}=[root,fdir(fi).name];
    end
end
ea_seed2map_weighted('/autofs/cluster/nimlab/USERS/Andy/MNI_connectomes/HCP_MGH_30fold_groupconnectome_gqi_lite.mat',...
    fis,...
    '/autofs/cluster/nimlab/USERS/Andy/MNI_connectomes/template_to_use.nii',maxdist);


exit
