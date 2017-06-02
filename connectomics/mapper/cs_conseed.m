function cs_conseed(dofMRI,dodMRI,dfold,sfile,cmd,writeoutsinglefiles,outputfolder,outputmask,dmrispace,fcname,dcname)
% wrapper for both dmri and fmri to generate seed2map files

if strcmp(outputmask,'.')
   outputmask=[]; 
end
if ~exist('dmrispace','var')
   dmrispace='222'; 
end
if ~isdeployed
    addpath(genpath('/autofs/cluster/nimlab/connectomes/software/lead_dbs'));
    addpath('/autofs/cluster/nimlab/connectomes/software/spm12');
    addpath(genpath('/home/agh14/lead_dbs'));
    addpath('/home/agh14/spm12');
end

fcdefault='GSP 1000 (Yeo 2011)';
dcdefault='HCP_MGH_30fold_groupconnectome (Horn 2017)';

if ~exist('fcname','var')
    fcname=fcdefault;
else
    if isempty(fcname)
        fcname=fcdefault;
    end
end

fcname=strrep(fcname,':','>');

if ~exist('dcname','var')
    dcname=dcdefault;
else
    if isempty(dcname)
        dcname=dcdefault;
    end
end
dcname=strrep(dcname,':','>');

if ~strcmp(outputfolder(end),filesep)
   outputfolder=[outputfolder,filesep]; 
end

% reformat inputs ? everything will be received as a string.
if ischar(dofMRI)
    dofMRI=str2double(dofMRI);
end
if ischar(dodMRI)
    dodMRI=str2double(dodMRI);
end
if ischar(writeoutsinglefiles)
    writeoutsinglefiles=str2double(writeoutsinglefiles);
end

if dofMRI
    ndfold=[dfold,filesep];
    cs_fmri_conseed(dfold,fcname,sfile,cmd,writeoutsinglefiles,outputfolder,outputmask);
end

if dodMRI
    ndfold=[dfold,filesep];
    cs_dmri_conseed(dfold,dcname,sfile,cmd,writeoutsinglefiles,outputfolder,outputmask,dmrispace);
end

exit