function ea_generate_datasetDescription(dest_filepath,flag,postop_modality)
[parent_dir,parent_file] = fileparts(dest_filepath);
dataset_description.Name = parent_file;
dataset_description.BIDSVersion = '1.6.0';
dataset_description.LEADVersion = '2.6';
dataset_description.DatasetType = 'raw'; %for backwards compatibility, as suggested by BIDS (by default)

if exist('flag', 'var') && strcmp(flag, 'raw')
   dataset_description.DatasetType = 'rawdata';
   dataset_description.Sessions = {'ses-preop','ses-postop'};
   dataset_description.preop_modality = 'MRI';
   if exist('postop_modality','var')
       dataset_description.postop_modality = postop_modality;
       dataset_description.HowToAcknowledge = '';
   end
end

output_file = fullfile(parent_dir,parent_file,'dataset_description.json');
json_fid = fopen(output_file,'w');
encodejson = savejson('',dataset_description);
fprintf(json_fid,encodejson);

if strcmp(flag,'root_folder')
    generate_bidsIgnore(parent_dir,parent_file);
end

%     if exist('flag','var') && strcmp(flag,'derivatives')
%         [filepath,filename,ext] = fileparts(dest_filepath);
%         description = generate_pipelineDescription(filename);
%         dataset_description.DatasetType = 'derivatives'; 
%         dataset_description.GeneratedBy.Name = filename;
%         dataset_description.GeneratedBy.Description = description{1};
%         dataset_description.HowToAcknowledge = '';
%     end
    
    
function description = generate_pipelineDescription(which_pipeline)
if strcmp(which_pipeline,'coregistration')
    description = {'Co-registration performs linear registration between the anchor modality and the anat modality of your choice'};
elseif strcmp(which_pipeline,'normalization')
    description = {'Normalization performs registration between the anchor modality(pre & post-op) and the standard space (MNI152NLin2009bAsym space) using registration algorithm specified'};
elseif strcmp(which_pipeline,'brainshift')
    description = {'Brainshift correction refines the fit of subcortical structures by correcting between preop T1w and postop CT'};
elseif strcmp(which_pipeline,'reconstruction')
    description = {'Reconstruction contains information about the localization of electrodes in the postop modality file'};
elseif strcmp(which_pipeline,'prefs')
    description = {'This folder contains user prefrences'};
elseif strcmp(which_pipeline,'log')
    description = {'This folder contains subject wide logs such as entire pipeline methods'};
elseif strcmp(which_pipeline,'miscellaneous')
    description = {'This folder contains files generated from earlier versions of lead-DBS, no longer required by the main pipeline'};
elseif strcmp(which_pipeline,'preprocessing')
    description = {'This folder contains minimally processed raw data - cropped and resliced used for running lead-DBS'};
elseif strcmp(which_pipeline,'stimulations')
    description = {'This folder contains Volume of tissue activated (vta) and E-fields files in MNI152NLin2009bAsym space as well as in native space'};
elseif strcmp(which_pipeline,'headmodel')
    description = {'This folder contains headmodel.mat files'};
else
    description = {''};
end


function generate_bidsIgnore(parent_dir,parent_file)
output_file_ignore = fullfile(parent_dir,parent_file,'.bidsignore');
opfile_ignore_fid = fopen(output_file_ignore,'w');
ignore_files = {'*CT*','*secondstepmask.nii','*thirdstepmask','*brainmask','*anchorNative*','*MNI152NLin2009bAsym*',...
    '*preproc*','*brainshiftmethod.json','*.xlsx','coregMR*','coregCT*','normalize*','*method.json','*ea_*.txt','*.mat',...
    '*rawimages*','export','miscellaneous','stimulations','atlases'};
for files = 1:length(ignore_files)
    ignore_filenames = [ignore_files{files} '\n'];
    fprintf(opfile_ignore_fid,ignore_filenames);
end
