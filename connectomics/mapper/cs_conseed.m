function cs_conseed(dofMRI,dodMRI,dfold,sfile,cmd,writeoutsinglefiles,outputfolder,outputmask)
% wrapper for both dmri and fmri to generate seed2map files

if ~isdeployed
    addpath(genpath('/autofs/cluster/nimlab/connectomes/software/lead_dbs'));
    addpath('/autofs/cluster/nimlab/connectomes/software/spm12');
end

if ~strcmp(outputfolder(end),filesep)
   outputfolder=[outputfolder,filesep]; 
end

if dofMRI
    ndfold=[dfold,filesep];
    cs_fmri_conseed(dfold,sfile,cmd,writeoutsinglefiles,outputfolder,outputmask);
end

if dodMRI
    ndfold=[dfold,filesep];
    cs_dmri_conseed([dfold,'dMRI',filesep,'HCP_MGH_30fold_groupconnectome_gqi_lite.mat'],sfile,cmd,writeoutsinglefiles,outputfolder,outputmask);
end