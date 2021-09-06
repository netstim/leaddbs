function results = legacy2bids(source,dest,isdicom,dicom_source)
%SANTIY: first check the existence of your source and dest directory. If it is not there,
%create a new directory.
if ~exist(source,'dir')
    warndlg("The source directory you have specified does not exist")
else
    addpath(source)
end
if ~exist(dest, 'dir')
    addpath(dest)
    mkdir(dest)
else
    addpath(dest)
end
%if you have dicom files and have not provided a dicom source directory
%(this can also just be source directory,so refactor this) then throw
%an error
if isdicom && ~exist(dicom_source,'var')
    warndlg("You have specified that you want dicom import, but have not specified a dicom source file. Please specific your dicom source directory!")
end

%define names of the new directorey structure
modes = {'anat','func'};
sessions = {'ses-preop','ses-postop'};
dicom_sessions = {'preop_mri','postop_ct'};
subfolder_cell = {'sourcedata','derivatives','rawdata'};
pipelines = {'brainshift','coregistration','normalization','reconstruction','preprocessing','prefs','log','export'};
%not sure how to handle log and export yet


%[[[currently creating a script that will take an entire source directory
%such as: dataset/{sub_01,sub_02,sub_03,..}
%get all subject names. This part of the code will be restructured when
%we get the dataset fetcher.]]]

%dir without dots is a modified dir command that with return the
%directory contents without the dots. those evaluate to folders and
%cause issues.
%by default, use the sub- prefix. If the patient already has a sub-, this
%flag is set to zero. 
flag = 1;
files = dir_without_dots(source);
sub_names = {files.name};
%get all patient names, again this will be changed based on the dataset
%fetcher command
for subs = 1:length(sub_names)
    patient_indx = 1;
    if isfolder(fullfile(source,sub_names{subs}))
        
        patient_names{patient_indx} = sub_names{subs};
        if startsWith(patient_names{patient_indx},'sub-')
            flag = 0;
        end
        patient_indx = patient_indx+1;
    end
end
%mapping will allow quick reference of the files to move
[brainshift_corr,coreg,normalization,preproc,recon,prefs] = create_bids_mapping();
%perform moving for each pat
for patients = 1:length(patient_names)
    
    patient = patient_names{patients};
    files_in_pat_folder = dir_without_dots(fullfile(source,patient));
    file_names = {files_in_pat_folder.name};
    file_index = 1;
    for i=1:length(file_names)
        %if they are already files, then add their full path to the
        %filelist to be moved
        if endsWith(file_names{i},'.nii') || endsWith(file_names{i},'.mat') || endsWith(file_names{i},'.h5') || endsWith(file_names{i},'.txt') || endsWith(file_names{i},'.log') || endsWith(file_names{i},'.gz') 
            files_to_move{file_index} = file_names{i};
            file_index = file_index + 1;
        end
        %now let's deal with subfolders
    end
    dir_names = {files_in_pat_folder([files_in_pat_folder.isdir]).name};
    for i=1:length(dir_names)
        %%%leave atlases as it is, it will be moved directly
        if ~(strcmp(dir_names{i},'atlases') || strcmp(dir_names{i},'stimulations'))
            this_folder = dir_without_dots(fullfile(source,patient,dir_names{i}));
            file_in_this_folder= {this_folder.name};
            file_inside_dir_indx = 1;
            for k=1:length(file_in_this_folder)
                if endsWith(file_in_this_folder{k},'.nii') || endsWith(file_in_this_folder{k},'.mat') || endsWith(file_in_this_folder{k},'.h5') || endsWith(file_in_this_folder{k},'.txt') || endsWith(file_in_this_folder{k},'.png')
                    if exist('files_to_move','var')
                        files_to_move{end+1} = file_in_this_folder{k};
                    else
                        files_to_move{file_inside_dir_indx} = file_in_this_folder{k};
                        file_inside_dir_indx = file_inside_dir_indx + 1; 
                    end
                end
            end
       
        else
            movefile(fullfile(source,patient,dir_names{i}),fullfile(dest,'derivatives','leaddbs',patient,dir_names{i}))
        end
    end
    %so we have created a list of the files to move, now we can create
    %the new dirrectory structure and move the correct files into the
    %correct "BIDS" directory
    for subfolders = 1:length(subfolder_cell)
        switch subfolder_cell{subfolders}
            case 'rawdata'
                for i=1:length(sessions)
                    new_path = fullfile(dest,subfolder_cell{subfolders},patient,sessions{i});
                    mkdir(new_path)
                    if strcmp(modes{i},'anat') && strcmp(sessions{i},'ses-preop')
                        %files to be moved into pre-op:raw_anat_*.nii
                        for files=1:length(files_to_move)
                            if regexp(files_to_move{files},'raw_anat_.*.nii')
                                movefile(fullfile(source,patient,files_to_move{files}),fullfile(new_path,modes{folders}))
                            end
                        end
                    elseif strcmp(modes{i},'anat') && strcmp(sessions{i},'ses-postop')
                        for files=1:length(files_to_move)
                            if ~isempty(regexp(files_to_move{files},'raw.*postop.*.nii','once')) || ~strcmp(files_to_move{files},'postop_ct.nii')
                                movefile(fullfile(source,patient,files_to_move{files}),fullfile(new_path,modes{folders}))
                            end
                        end
                      
                    end
                end
            case 'derivatives'
                new_path = fullfile(dest,subfolder_cell{subfolders},'leaddbs',patient);
                for folders=1:length(pipelines)
                    mkdir(fullfile(new_path,pipelines{folders}))
                    for files=1:length(files_to_move)
                        if strcmp(pipelines{folders},'coregistration')
                            mkdir(fullfile(new_path,pipelines{folders},'anat'));
                            mkdir(fullfile(new_path,pipelines{folders},'log'));
                            mkdir(fullfile(new_path,pipelines{folders},'checkreg'));
                            mkdir(fullfile(new_path,pipelines{folders},'transformations'))
                            if ismember(files_to_move{files},coreg{:,1})
                                %if ~isempty(regexp(files_to_move{files},'^anat_t[1,2].nii||^postop_.*nii||.*rpostop.*||.*coreg.*') == 1)
                                %corresponding index of the new pat
                                indx = cellfun(@(x)strcmp(x,files_to_move{files}),coreg{:,1});
                                if endsWith(files_to_move{files},'.nii')
                                    %first move%
                                    movefile(fullfile(source,patient,files_to_move{files}),fullfile(new_path,pipelines{folders},'anat',files_to_move{files}));
                                    %then rename%
                                    if flag
                                        movefile(fullfile(new_path,pipelines{folders},'anat',files_to_move{files}),fullfile(new_path,pipelines{folders},'anat',['sub-','_',patient,coreg{1,2}{indx}]));
                                    else
                                         movefile(fullfile(new_path,pipelines{folders},'anat',files_to_move{files}),fullfile(new_path,pipelines{folders},'anat',[patient,'_',coreg{1,2}{indx}]));
                                    end
                                elseif endsWith(files_to_move{files},'.png')
                                    movefile(fullfile(source,patient,'checkreg',files_to_move{files}),fullfile(new_path,pipelines{folders},'checkreg',files_to_move{files}));
                                    if flag
                                        movefile(fullfile(new_path,pipelines{folders},'checkreg',files_to_move{files}),fullfile(new_path,pipelines{folders},'checkreg',['sub-','_',patient,coreg{1,2}{indx}]));
                                    else
                                         movefile(fullfile(new_path,pipelines{folders},'checkreg',files_to_move{files}),fullfile(new_path,pipelines{folders},'checkreg',[patient,'_',coreg{1,2}{indx}]));
                                    end
                                elseif endsWith('.log') || ~isempty(regexp(files_to_move{files},'.*_approved||.*_applied.mat'))
                                    movefile(fullfile(source,patient,files_to_move{files}),fullfile(new_path,pipelines{folders},'log',files_to_move{files}));
                                    if flag
                                        movefile(fullfile(new_path,pipelines{folders},'log',files_to_move{files}),fullfile(new_path,pipelines{folders},'log',['sub-','_',patient,coreg{1,2}{indx}]));
                                    else
                                         movefile(fullfile(new_path,pipelines{folders},'log',files_to_move{files}),fullfile(new_path,pipelines{folders},'log',[patient,'_',coreg{1,2}{indx}]));
                                    end
                                elseif endsWith(files_to_move{files},'.mat')
                                   movefile(fullfile(source,patient,files_to_move{files}),fullfile(new_path,pipelines{folders},'transformations',files_to_move{files}));
                                   if flag
                                        movefile(fullfile(new_path,pipelines{folders},'transformations',files_to_move{files}),fullfile(new_path,pipelines{folders},'transformations',['sub-','_',patient,coreg{1,2}{indx}]));
                                    else
                                         movefile(fullfile(new_path,pipelines{folders},'transformations',files_to_move{files}),fullfile(new_path,pipelines{folders},'transformations',[patient,'_',coreg{1,2}{indx}]));
                                    end
                                end
                            elseif ~isempty(regexp(files_to_move{files},'^coreg.*.log')) 
                                    movefile(fullfile(source,patient,files_to_move{files}),fullfile(new_path,pipelines{folders},'log',files_to_move{files}));
                                    
                            end
         
                        elseif strcmp(pipelines{folders},'brainshift')
                            mkdir(fullfile(new_path,pipelines{folders},'anat'));
                            mkdir(fullfile(new_path,pipelines{folders},'log'));
                            mkdir(fullfile(new_path,pipelines{folders},'checkreg'));
                            mkdir(fullfile(new_path,pipelines{folders},'transformations'))
                            if ismember(files_to_move{files},brainshift_corr{:,1})
                                indx = cellfun(@(x)strcmp(x,files_to_move{files}),brainshift_corr{:,1});
                                if endsWith(files_to_move{files},'.nii')
                                    %first move%
                                    movefile(fullfile(source,patient,'scrf',files_to_move{files}),fullfile(new_path,pipelines{folders},'anat',files_to_move{files}));
                                    %then rename%
                                    if flag
                                        movefile(fullfile(new_path,pipelines{folders},'anat',files_to_move{files}),fullfile(new_path,pipelines{folders},'anat',['sub-','_',patient,brainshift_corr{1,2}{indx}]));
                                    else
                                        movefile(fullfile(new_path,pipelines{folders},'anat',files_to_move{files}),fullfile(new_path,pipelines{folders},'anat',[patient,'_',brainshift_corr{1,2}{indx}]));
                                    end
                                elseif endsWith(files_to_move{files},'.png')
                                    movefile(fullfile(source,patient,'scrf',files_to_move{files}),fullfile(new_path,pipelines{folders},'checkreg',files_to_move{files}));
                                    if flag
                                        movefile(fullfile(new_path,pipelines{folders},'checkreg',files_to_move{files}),fullfile(new_path,pipelines{folders},'checkreg',['sub-','_',patient,brainshift_corr{1,2}{indx}]));
                                    else
                                         movefile(fullfile(new_path,pipelines{folders},'checkreg',files_to_move{files}),fullfile(new_path,pipelines{folders},'checkreg',[patient,'_',brainshift_corr{1,2}{indx}]));
                                    end
                                elseif endsWith(files_to_move{files},'.txt') || ~isempty(regexp(files_to_move{files},'.*_approved||.*_applied.mat'))
                                    movefile(fullfile(source,patient,'scrf',files_to_move{files}),fullfile(new_path,pipelines{folders},'log',files_to_move{files}));
                                   if flag
                                        movefile(fullfile(new_path,pipelines{folders},'log',files_to_move{files}),fullfile(new_path,pipelines{folders},'log',['sub-','_',patient,brainshift_corr{1,2}{indx}]));
                                    else
                                         movefile(fullfile(new_path,pipelines{folders},'log',files_to_move{files}),fullfile(new_path,pipelines{folders},'log',[patient,'_',brainshift_corr{1,2}{indx}]));
                                    end
                                elseif endsWith(files_to_move{files},'.mat')
                                   movefile(fullfile(source,patient,'scrf',files_to_move{files}),fullfile(new_path,pipelines{folders},'transformations',files_to_move{files}));
                                   if flag
                                        movefile(fullfile(new_path,pipelines{folders},'transformations',files_to_move{files}),fullfile(new_path,pipelines{folders},'transformations',['sub-','_',patient,brainshift_corr{1,2}{indx}]));
                                    else
                                         movefile(fullfile(new_path,pipelines{folders},'transformations',files_to_move{files}),fullfile(new_path,pipelines{folders},'transformations',[patient,'_',brainshift_corr{1,2}{indx}]));
                                    end
                                end
                            end
                        elseif strcmp(pipelines{folders},'normalization')
                            mkdir(fullfile(new_path,pipelines{folders},'anat'));
                            mkdir(fullfile(new_path,pipelines{folders},'log'));
                            mkdir(fullfile(new_path,pipelines{folders},'checkreg'));
                            mkdir(fullfile(new_path,pipelines{folders},'transformations'))
                            if ismember(files_to_move{files},normalization{:,1})
                                %corresponding index of the new pat
                                indx = cellfun(@(x)strcmp(x,files_to_move{files}),normalization{:,1});
                                if endsWith(files_to_move{files},'.nii')
                                    %first move%
                                    movefile(fullfile(source,patient,files_to_move{files}),fullfile(new_path,pipelines{folders},'anat',files_to_move{files}));
                                    %then rename%
                                    if flag
                                        movefile(fullfile(new_path,pipelines{folders},'anat',files_to_move{files}),fullfile(new_path,pipelines{folders},'anat',['sub-','_',patient,normalization{1,2}{indx}]));
                                    else
                                        movefile(fullfile(new_path,pipelines{folders},'anat',files_to_move{files}),fullfile(new_path,pipelines{folders},'anat',[patient,'_',normalization{1,2}{indx}]));
                                    end
                                elseif endsWith(files_to_move{files},'.png')
                                    movefile(fullfile(source,patient,'checkreg',files_to_move{files}),fullfile(new_path,pipelines{folders},'checkreg',files_to_move{files}));
                                    if flag
                                        movefile(fullfile(new_path,pipelines{folders},'checkreg',files_to_move{files}),fullfile(new_path,pipelines{folders},'checkreg',['sub-','_',patient,normalization{1,2}{indx}]));
                                    else
                                        movefile(fullfile(new_path,pipelines{folders},'checkreg',files_to_move{files}),fullfile(new_path,pipelines{folders},'checkreg',[patient,'_',normalization{1,2}{indx}]));
                                    end
                                elseif endsWith(files_to_move{files},'.log') || ~isempty(regexp(files_to_move{files},'.*_approved||.*_applied.mat'))
                                   movefile(fullfile(source,patient,files_to_move{files}),fullfile(new_path,pipelines{folders},'log',files_to_move{files}));
                                   if flag
                                        movefile(fullfile(new_path,pipelines{folders},'log',files_to_move{files}),fullfile(new_path,pipelines{folders},'log',['sub-','_',patient,normalization{1,2}{indx}]));
                                    else
                                         movefile(fullfile(new_path,pipelines{folders},'log',files_to_move{files}),fullfile(new_path,pipelines{folders},'log',[patient,'_',normalization{1,2}{indx}]));
                                    end
                                elseif endsWith(files_to_move{files},'.mat') || endsWith(files_to_move{files},'.gz')
                                   movefile(fullfile(source,patient,files_to_move{files}),fullfile(new_path,pipelines{folders},'transformations',files_to_move{files}));
                                   if flag
                                        movefile(fullfile(new_path,pipelines{folders},'transformations',files_to_move{files}),fullfile(new_path,pipelines{folders},'transformations',['sub-','_',patient,normalization{1,2}{indx}]));
                                    else
                                        movefile(fullfile(new_path,pipelines{folders},'transformations',files_to_move{files}),fullfile(new_path,pipelines{folders},'transformations',[patient,'_',normalization{1,2}{indx}]));
                                    end
                                end
                            elseif ~isempty(regexp(files_to_move{files},'^normalize_.*.log'))
                                movefile(fullfile(source,patient,files_to_move{files}),fullfile(new_path,pipelines{folders},'log',files_to_move{files}));
                                
                            end
                        elseif strcmp(pipelines{folders},'reconstruction')
                            mkdir(fullfile(new_path,pipelines{folders}));
                            if ismember(files_to_move{files},recon{:,1})
                                %corresponding index of the new pat
                                indx = cellfun(@(x)strcmp(x,files_to_move{files}),recon{:,1});
                                movefile(fullfile(source,patient,files_to_move{files}),fullfile(new_path,pipelines{folders},files_to_move{files}));
                                if flag
                                    movefile(fullfile(new_path,pipelines{folders},files_to_move{files}),fullfile(new_path,pipelines{folders},['sub-','_',patient,recon{1,2}{indx}]));
                                else
                                    movefile(fullfile(new_path,pipelines{folders},files_to_move{files}),fullfile(new_path,pipelines{folders},[patient,'_',recon{1,2}{indx}]));
                                end
                                
                            end
                        elseif strcmp(pipelines{folders},'preprocessing')
                            mkdir(fullfile(new_path,pipelines{folders},'anat'));
                            mkdir(fullfile(new_path,pipelines{folders},'log'));
                            if ismember(files_to_move{files},preproc{:,1})
                                %corresponding index of the new pat
                                indx = cellfun(@(x)strcmp(x,files_to_move{files}),preproc{:,1});
                                if endsWith(files_to_move{files},'.nii')
                                    %first move%
                                    movefile(fullfile(source,patient,files_to_move{files}),fullfile(new_path,pipelines{folders},'anat',files_to_move{files}));
                                    %then rename%
                                    if flag
                                        movefile(fullfile(new_path,pipelines{folders},'anat',files_to_move{files}),fullfile(new_path,pipelines{folders},'anat',['sub-','_',patient,preproc{1,2}{indx}]));
                                    else
                                        movefile(fullfile(new_path,pipelines{folders},'anat',files_to_move{files}),fullfile(new_path,pipelines{folders},'anat',[patient,'_',preproc{1,2}{indx}]));
                                    end
                                elseif endsWith(files_to_move{files},'.log')
                                    movefile(fullfile(source,patient,files_to_move{files}),fullfile(new_path,pipelines{folders},'log',files_to_move{files}));
                                    if flag
                                        movefile(fullfile(new_path,pipelines{folders},'log',files_to_move{files}),fullfile(new_path,pipelines{folders},'log',['sub-','_',patient,preproc{1,2}{indx}]));
                                    else
                                         movefile(fullfile(new_path,pipelines{folders},'log',files_to_move{files}),fullfile(new_path,pipelines{folders},'log',[patient,'_',preproc{1,2}{indx}]));
                                    end
                                end
                            end
                        elseif strcmp(pipelines{folders},'prefs')
                            if ismember(files_to_move{files},prefs{:,1})
                                disp(files_to_move{files})
                                %corresponding index of the new pat
                                indx = cellfun(@(x)strcmp(x,files_to_move{files}),prefs{:,1});
                                movefile(fullfile(source,patient,files_to_move{files}),fullfile(new_path,pipelines{folders},files_to_move{files}));
                                if flag
                                    movefile(fullfile(new_path,pipelines{folders},files_to_move{files}),fullfile(new_path,pipelines{folders},['sub-','_',patient,prefs{1,2}{indx}]));
                                else
                                    movefile(fullfile(new_path,pipelines{folders},files_to_move{files}),fullfile(new_path,pipelines{folders},[patient,'_',prefs{1,2}{indx}]));
                                end
                            end
                         end
                     end
                 end
            otherwise
                for i=1:length(dicom_sessions)
                    new_path = fullfile(dest,subfolder_cell{subfolders},patient,'DICOM',dicom_sessions{i});
                    if ~exist(new_path,'dir')
                        mkdir(new_path)
                    end
                    if ~isdicom
                        disp("There are no dicom images, source data folder will be empty")
                    else
                        movefile(dicom_source,new_path)
                    end
                end
         end
     end
 end
%     %writing another script which will return the mapping between the
%     %legacy files and the new files.
%    % [old_files,new_files] = create_bids_mapping;
%     %comparing the files in our list of files and renaming the files into
%     %the correct new name
end

