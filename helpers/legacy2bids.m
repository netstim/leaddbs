function results = legacy2bids(source,dest,isdicom,dicom_source,move)
%SANTIY: first check the existence of your source and dest directory. If it is not there,
%create a new directory.
for j=1:length(source)
    if ~exist(source{j},'dir')
        warndlg("The source directory you have specified does not exist")
    else
        addpath(source{j})
    end
    if ~exist(dest, 'var')
        dest{j} = source{j};
        addpath(dest{j});
    elseif exist(dest{j},'var') && ~exist(dest{j},'dir')
        addpath(dest{j});
        mkdir(dest{j});
    end
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
subfolder_cell = {'sourcedata','rawdata','derivatives'};
pipelines = {'brainshift','coregistration','normalization','reconstruction','preprocessing','prefs','log','export'};
%not sure how to handle log and export yet
%mapping will allow quick reference of the files to move
[brainshift_corr,coreg,normalization,preproc,recon,prefs] = create_bids_mapping();
%dir without dots is a modified dir command that with return the
%directory contents without the dots. those evaluate to folders and
%cause issues.
%by default, use the sub- prefix. If the patient already has a sub-, this
%flag is set to zero.

%get all patient names, again this will be changed based on the dataset
%fetcher command
%perform moving for each pat
for patients = 1:length(source)
    source_patient = source{patients};
    dest_patient = dest{patients};
    flag = 1;
    if startsWith(source_patient,'sub-')
        flag = 0;
        
    end
    [filepath,patient_name,ext] = fileparts(source_patient);
    
    
    files_in_pat_folder = dir_without_dots(source_patient);
    file_names = {files_in_pat_folder.name};
    file_index = 1;
    for j=1:length(file_names)
        %if they are already files, then add their full path to the
        %filelist to be moved
        if endsWith(file_names{j},'.nii') || endsWith(file_names{j},'.mat') || endsWith(file_names{j},'.h5') || endsWith(file_names{j},'.txt') || endsWith(file_names{j},'.log') || endsWith(file_names{j},'.gz')
            files_to_move{file_index} = file_names{j};
            file_index = file_index + 1;
        end
        %now let's deal with subfolders
    end
    dir_names = {files_in_pat_folder([files_in_pat_folder.isdir]).name};
    for j=1:length(dir_names)
        %%%leave atlases as it is, it will be moved directly
        if ~(strcmp(dir_names{j},'atlases') || strcmp(dir_names{j},'stimulations'))
            this_folder = dir_without_dots(fullfile(source_patient,dir_names{j}));
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
            if move
                movefile(fullfile(source_patient,dir_names{j}),fullfile(dest_patient,'derivatives','leaddbs',patient_name,dir_names{j}));
            else
                copyfile(fullfile(source_patient,dir_names{j}),fullfile(dest_patient,'derivatives','leaddbs',patient_name));
            end
        end
    end
    %so we have created a list of the files to move, now we can create
    %the new dirrectory structure and move the correct files into the
    %correct "BIDS" directory
    for subfolders = 1:length(subfolder_cell)
        switch subfolder_cell{subfolders}
            case 'rawdata'
                for i=1:length(modes)
                    for j=1:length(sessions)
                        new_path = fullfile(dest_patient,subfolder_cell{subfolders},patient_name,sessions{j});
                        if ~exist(new_path,'dir')
                            mkdir(new_path)
                        end
                        for files=1:length(files_to_move)
                            if strcmp(modes{i},'anat') && strcmp(sessions{j},'ses-preop')
                                %files to be moved into pre-op:raw_anat_*.nii
                                if regexp(files_to_move{files},'raw_anat_.*.nii')
                                    if exist(fullfile(source_patient,files_to_move{files}),'file')
                                        str_var = strsplit(files_to_move{files},'_');
                                        modality_str = strsplit(str_var{end},'.');
                                        if strcmp(modality_str{1},'t1')
                                            modality = 'T1w';
                                        elseif strcmp(modality_str{1},'t2')
                                            modality = 'T2w';
                                        else
                                            modality = upper(modality_str{1});
                                        end
                                        if move
                                            movefile(fullfile(source_patient,files_to_move{files}),fullfile(new_path,files_to_move{files}));
                                        else
                                            copyfile(fullfile(source_patient,files_to_move{files}),new_path);
                                        end
                                        if flag
                                            movefile(fullfile(new_path,files_to_move{files}),fullfile(new_path,['sub-',patient_name,'_',sessions{j},'_',modality,'.nii']));
                                        else
                                            movefile(fullfile(new_path,files_to_move{files}),fullfile(new_path,[patient_name,'_',modality,'.nii']));
                                        end
                                    end
                                end
                            elseif strcmp(modes{i},'anat') && strcmp(sessions{j},'ses-postop')
                                if ~isempty(regexp(files_to_move{files},'raw_postop_.*.nii')) || strcmp(files_to_move{files},'postop_ct.nii')
                                    if exist(fullfile(source_patient,files_to_move{files}),'file')
                                        modality_str = strsplit(files_to_move{files},'_');
                                        if strcmp(modality_str{end},'t1')
                                            modality = 'T1w';
                                        elseif strcmp(modality_str{end},'t2')
                                            modality = 'T2w';
                                        elseif strcmp(modality_str{end},'CT')
                                            modality = 'CT';
                                        else
                                            str_var = strsplit(modality_str{end},'.');
                                            modality = upper(str_var{1});
                                        end
                                        if move
                                            movefile(fullfile(source_patient,files_to_move{files}),fullfile(new_path,files_to_move{files}));
                                        else
                                            copyfile(fullfile(source_patient,files_to_move{files}),new_path);
                                        end
                                        if flag
                                            movefile(fullfile(new_path,files_to_move{files}),fullfile(new_path,['sub-',patient_name,'_',sessions{j},'_',modality,'.nii']));
                                        else
                                            movefile(fullfile(new_path,files_to_move{files}),fullfile(new_path,[patient_name,'_',modality,'.nii']));
                                        end
                                    end
                                end
                            end
                            
                        end
                    end
                end
            case 'derivatives'
                new_path = fullfile(dest_patient,subfolder_cell{subfolders},'leaddbs',patient_name);
                for folders=1:length(pipelines)
                    if ~exist(fullfile(new_path,pipelines{folders}),'dir')
                        mkdir(fullfile(new_path,pipelines{folders}))
                    end
                    for files=1:length(files_to_move)
                        if strcmp(pipelines{folders},'coregistration')
                            anat_dir = fullfile(new_path,pipelines{folders},'anat');
                            log_dir = fullfile(new_path,pipelines{folders},'log');
                            checkreg_dir = fullfile(new_path,pipelines{folders},'checkreg');
                            transformations_dir = fullfile(new_path,pipelines{folders},'transformations');
                            if ~exist(anat_dir,'dir')
                                mkdir(anat_dir)
                            end
                            if ~exist(log_dir,'dir')
                                mkdir(log_dir)
                            end
                            if ~exist(checkreg_dir,'dir')
                                mkdir(checkreg_dir)
                            end
                            if ~exist(transformations_dir,'dir')
                                mkdir(transformations_dir)
                            end
                            if ismember(files_to_move{files},coreg{:,1})
                                %if ~isempty(regexp(files_to_move{files},'^anat_t[1,2].nii||^postop_.*nii||.*rpostop.*||.*coreg.*') == 1)
                                %corresponding index of the new pat
                                indx = cellfun(@(x)strcmp(x,files_to_move{files}),coreg{:,1});
                                if endsWith(files_to_move{files},'.nii')
                                    %first move%
                                    if exist(fullfile(source_patient,files_to_move{files}),'file')
                                        if move
                                            movefile(fullfile(source_patient,files_to_move{files}),fullfile(anat_dir,files_to_move{files}));
                                        else
                                           copyfile(fullfile(source_patient,files_to_move{files}),anat_dir);  
                                        end
                                        if flag
                                            movefile(fullfile(new_path,pipelines{folders},'anat',files_to_move{files}),fullfile(anat_dir,['sub-','_',patient_name,coreg{1,2}{indx}]));
                                        else
                                            movefile(fullfile(new_path,pipelines{folders},'anat',files_to_move{files}),fullfile(anat_dir,[patient_name,'_',coreg{1,2}{indx}]));
                                        end
                                    elseif exist(fullfile(source_patient,'scrf',files_to_move{files}),'file')
                                        if move
                                            movefile(fullfile(source_patient,'scrf',files_to_move{files}),fullfile(anat_dir,files_to_move{files}));
                                        else
                                            copyfile(fullfile(source_patient,'scrf',files_to_move{files}),anat_dir)
                                        end
                                        %then rename%
                                        if flag
                                            movefile(fullfile(new_path,pipelines{folders},'anat',files_to_move{files}),fullfile(anat_dir,['sub-','_',patient_name,coreg{1,2}{indx}]));
                                        else
                                            movefile(fullfile(new_path,pipelines{folders},'anat',files_to_move{files}),fullfile(anat_dir,[patient_name,'_',coreg{1,2}{indx}]));
                                        end
                                    end
                                elseif endsWith(files_to_move{files},'.png')
                                    if exist(fullfile(source_patient,'checkreg',files_to_move{files}),'file')
                                        if move
                                            movefile(fullfile(source_patient,'checkreg',files_to_move{files}),fullfile(checkreg_dir,files_to_move{files}));
                                        else
                                            copyfile(fullfile(source_patient,'checkreg',files_to_move{files}),checkreg_dir);
                                        end
                                        if flag
                                            movefile(fullfile(new_path,pipelines{folders},'checkreg',files_to_move{files}),fullfile(checkreg_dir,['sub-','_',patient_name,coreg{1,2}{indx}]));
                                        else
                                            movefile(fullfile(new_path,pipelines{folders},'checkreg',files_to_move{files}),fullfile(checkreg_dir,[patient_name,'_',coreg{1,2}{indx}]));
                                        end
                                    end
                                elseif endsWith(files_to_move{files},'.log') || ~isempty(regexp(files_to_move{files},'.*_approved||.*_applied.mat'))
                                    if exist(fullfile(source_patient,files_to_move{files}),'file')
                                        if move
                                            movefile(fullfile(source_patient,files_to_move{files}),fullfile(log_dir,files_to_move{files}));
                                        else
                                            copyfile(fullfile(source_patient,files_to_move{files}),log_dir)
                                        end
                                        if flag
                                            movefile(fullfile(new_path,pipelines{folders},'log',files_to_move{files}),fullfile(log_dir,['sub-','_',patient_name,coreg{1,2}{indx}]));
                                        else
                                            movefile(fullfile(new_path,pipelines{folders},'log',files_to_move{files}),fullfile(log_dir,[patient_name,'_',coreg{1,2}{indx}]));
                                        end
                                    end
                                elseif endsWith(files_to_move{files},'.mat')
                                    if exist(fullfile(source_patient,files_to_move{files}),'file')
                                        if move
                                            movefile(fullfile(source_patient,files_to_move{files}),fullfile(transformations_dir,files_to_move{files}));
                                        else
                                            copyfile(fullfile(source_patient,files_to_move{files}),transformations_dir)
                                        end
                                        if flag
                                            movefile(fullfile(new_path,pipelines{folders},'transformations',files_to_move{files}),fullfile(transformations_dir,['sub-','_',patient_name,coreg{1,2}{indx}]));
                                        else
                                            movefile(fullfile(new_path,pipelines{folders},'transformations',files_to_move{files}),fullfile(transformations_dir,[patient_name,'_',coreg{1,2}{indx}]));
                                        end
                                    end
                                end
                            elseif ~isempty(regexp(files_to_move{files},'^coreg.*.log'))
                                if exist(fullfile(source_patient,files_to_move{files}),'file')
                                    if move
                                        movefile(fullfile(source_patient,files_to_move{files}),fullfile(log_dir,files_to_move{files}));
                                    else
                                        copyfile(fullfile(source_patient,files_to_move{files}),log_dir)
                                    end
                                end
                            end
                            
                        elseif strcmp(pipelines{folders},'brainshift')
                            anat_dir = fullfile(new_path,pipelines{folders},'anat');
                            log_dir = fullfile(new_path,pipelines{folders},'log');
                            checkreg_dir = fullfile(new_path,pipelines{folders},'checkreg');
                            transformations_dir = fullfile(new_path,pipelines{folders},'transformations');
                            if ~exist(anat_dir,'dir')
                                mkdir(anat_dir)
                            end
                            if ~exist(log_dir,'dir')
                                mkdir(log_dir)
                            end
                            if ~exist(checkreg_dir,'dir')
                                mkdir(checkreg_dir)
                            end
                            if ~exist(transformations_dir,'dir')
                                mkdir(transformations_dir)
                            end
                            if ismember(files_to_move{files},brainshift_corr{:,1})
                                indx = cellfun(@(x)strcmp(x,files_to_move{files}),brainshift_corr{:,1});
                                if endsWith(files_to_move{files},'.nii')
                                    %first move%
                                    if exist(fullfile(source_patient,'scrf',files_to_move{files}),'file')
                                        if move
                                            movefile(fullfile(source_patient,'scrf',files_to_move{files}),fullfile(anat_dir,files_to_move{files}));
                                        else
                                            copyfile(fullfile(source_patient,'scrf',files_to_move{files}),anat_dir);
                                        end
                                        %then rename%
                                        if flag
                                            movefile(fullfile(new_path,pipelines{folders},'anat',files_to_move{files}),fullfile(anat_dir,['sub-',patient_name,'_',brainshift_corr{1,2}{indx}]));
                                        else
                                            movefile(fullfile(new_path,pipelines{folders},'anat',files_to_move{files}),fullfile(anat_dir,[patient_name,'_',brainshift_corr{1,2}{indx}]));
                                        end
                                    end
                                elseif endsWith(files_to_move{files},'.png')
                                    if exist(fullfile(source_patient,'scrf',files_to_move{files}),'file')
                                        if move
                                            movefile(fullfile(source_patient,'scrf',files_to_move{files}),fullfile(checkreg_dir,files_to_move{files}));
                                        else
                                            copyfile(fullfile(source_patient,'scrf',files_to_move{files}),checkreg_dir);
                                        end
                                        if flag
                                            movefile(fullfile(new_path,pipelines{folders},'checkreg',files_to_move{files}),fullfile(checkreg_dir,['sub-',patient_name,'_',brainshift_corr{1,2}{indx}]));
                                        else
                                            movefile(fullfile(new_path,pipelines{folders},'checkreg',files_to_move{files}),fullfile(checkreg_dir,[patient_name,'_',brainshift_corr{1,2}{indx}]));
                                        end
                                    end
                                elseif endsWith(files_to_move{files},'.txt') || ~isempty(regexp(files_to_move{files},'.*_approved||.*_applied.mat'))
                                    if exist(fullfile(source_patient,'scrf',files_to_move{files}),'file')
                                        if move
                                            movefile(fullfile(source_patient,'scrf',files_to_move{files}),fullfile(log_dir,files_to_move{files}));
                                        else
                                            copyfile(fullfile(source_patient,'scrf',files_to_move{files}),log_dir)
                                        end
                                        if flag
                                            movefile(fullfile(new_path,pipelines{folders},'log',files_to_move{files}),fullfile(log_dir,['sub-',patient_name,'_',brainshift_corr{1,2}{indx}]));
                                        else
                                            movefile(fullfile(new_path,pipelines{folders},'log',files_to_move{files}),fullfile(log_dir,[patient_name,'_',brainshift_corr{1,2}{indx}]));
                                        end
                                    end
                                elseif endsWith(files_to_move{files},'.mat')
                                    if exist(fullfile(source_patient,'scrf',files_to_move{files}),'file')
                                        if move
                                            movefile(fullfile(source_patient,'scrf',files_to_move{files}),fullfile(transformations_dir,files_to_move{files}));
                                        else
                                            copyfile(fullfile(source_patient,'scrf',files_to_move{files}),transformations_dir);
                                        end
                                        if flag
                                            movefile(fullfile(new_path,pipelines{folders},'transformations',files_to_move{files}),fullfile(transformations_dir,['sub-',patient_name,'_',brainshift_corr{1,2}{indx}]));
                                        else
                                            movefile(fullfile(new_path,pipelines{folders},'transformations',files_to_move{files}),fullfile(transformations_dir,[patient_name,'_',brainshift_corr{1,2}{indx}]));
                                        end
                                    end
                                end
                            end
                        elseif strcmp(pipelines{folders},'normalization')
                            anat_dir = fullfile(new_path,pipelines{folders},'anat');
                            log_dir = fullfile(new_path,pipelines{folders},'log');
                            checkreg_dir = fullfile(new_path,pipelines{folders},'checkreg');
                            transformations_dir = fullfile(new_path,pipelines{folders},'transformations');
                            if ~exist(anat_dir,'dir')
                                mkdir(anat_dir)
                            end
                            if ~exist(log_dir,'dir')
                                mkdir(log_dir)
                            end
                            if ~exist(checkreg_dir,'dir')
                                mkdir(checkreg_dir)
                            end
                            if ~exist(transformations_dir,'dir')
                                mkdir(transformations_dir)
                            end
                            if ismember(files_to_move{files},normalization{:,1})
                                %corresponding index of the new pat
                                indx = cellfun(@(x)strcmp(x,files_to_move{files}),normalization{:,1});
                                if endsWith(files_to_move{files},'.nii')
                                    if exist(fullfile(source_patient,files_to_move{files}),'file')
                                        %first move%
                                        if move
                                            movefile(fullfile(source_patient,files_to_move{files}),fullfile(anat_dir,files_to_move{files}));
                                        else
                                            copyfile(fullfile(source_patient,files_to_move{files}),anat_dir);
                                        end
                                        %then rename%
                                        if flag
                                            movefile(fullfile(new_path,pipelines{folders},'anat',files_to_move{files}),fullfile(anat_dir,['sub-',patient_name,'_',normalization{1,2}{indx}]));
                                        else
                                            movefile(fullfile(new_path,pipelines{folders},'anat',files_to_move{files}),fullfile(anat_dir,[patient_name,'_',normalization{1,2}{indx}]));
                                        end
                                    end
                                elseif endsWith(files_to_move{files},'.png')
                                    if exist(fullfile(source_patient,'checkreg',files_to_move{files}),'file')
                                        if move
                                            movefile(fullfile(source_patient,'checkreg',files_to_move{files}),fullfile(checkreg_dir,files_to_move{files}));
                                        else
                                            copyfile(fullfile(source_patient,'checkreg',files_to_move{files}),checkreg_dir)
                                        end
                                        if flag
                                            movefile(fullfile(new_path,pipelines{folders},'checkreg',files_to_move{files}),fullfile(checkreg_dir,['sub-',patient_name,'_',normalization{1,2}{indx}]));
                                        else
                                            movefile(fullfile(new_path,pipelines{folders},'checkreg',files_to_move{files}),fullfile(checkreg_dir,[patient_name,'_',normalization{1,2}{indx}]));
                                        end
                                    end
                                elseif endsWith(files_to_move{files},'.log') || ~isempty(regexp(files_to_move{files},'.*_approved||.*_applied.mat'))
                                    if exist(fullfile(source_patient,files_to_move{files}),'file')
                                        if move
                                            movefile(fullfile(source_patient,files_to_move{files}),fullfile(log_dir,files_to_move{files}));
                                        else
                                            copyfile(fullfile(source_patient,files_to_move{files}),log_dir);
                                        end
                                        
                                        if flag
                                            movefile(fullfile(new_path,pipelines{folders},'log',files_to_move{files}),fullfile(log_dir,['sub-',patient_name,'_',normalization{1,2}{indx}]));
                                        else
                                            movefile(fullfile(new_path,pipelines{folders},'log',files_to_move{files}),fullfile(log_dir,[patient_name,'_',normalization{1,2}{indx}]));
                                        end
                                    end
                                elseif endsWith(files_to_move{files},'.mat') || endsWith(files_to_move{files},'.gz')
                                    if exist(fullfile(source_patient,files_to_move{files}),'file')
                                        if move
                                            movefile(fullfile(source_patient,files_to_move{files}),fullfile(transformations_dir,files_to_move{files}));
                                        else
                                            copyfile(fullfile(source_patient,files_to_move{files}),transformations_dir)
                                        end
                                        if flag
                                            movefile(fullfile(new_path,pipelines{folders},'transformations',files_to_move{files}),fullfile(transformations_dir,['sub-',patient_name,'_',normalization{1,2}{indx}]));
                                        else
                                            movefile(fullfile(new_path,pipelines{folders},'transformations',files_to_move{files}),fullfile(transformations_dir,[patient_name,'_',normalization{1,2}{indx}]));
                                        end
                                    end
                                elseif ~isempty(regexp(files_to_move{files},'^normalize_.*.log')) || strcmp(files_to_move{files},'ea_ants_command.txt')
                                    if exist(fullfile(source_patient,files_to_move{files}),'file')
                                        if move
                                            movefile(fullfile(source_patient,files_to_move{files}),fullfile(log_dir,files_to_move{files}));
                                        else
                                            copyfile(fullfile(source_patient,files_to_move{files}),log_dir)
                                        end
                                    end
                                end
                            end
                        elseif strcmp(pipelines{folders},'reconstruction')
                            recon_dir = fullfile(new_path,pipelines{folders});
                            if ~exist(recon_dir,'dir')
                                mkdir(recon_dir)
                            end
                            if ismember(files_to_move{files},recon{:,1})
                                %corresponding index of the new pat
                                indx = cellfun(@(x)strcmp(x,files_to_move{files}),recon{:,1});
                                if move
                                    movefile(fullfile(source_patient,files_to_move{files}),fullfile(recon_dir,files_to_move{files}));
                                else
                                    copyfile(fullfile(source_patient,files_to_move{files}),recon_dir)
                                end
                                
                                if flag
                                    movefile(fullfile(new_path,pipelines{folders},files_to_move{files}),fullfile(recon_dir,['sub-',patient_name,'_',recon{1,2}{indx}]));
                                else
                                    movefile(fullfile(new_path,pipelines{folders},files_to_move{files}),fullfile(recon_dir,[patient_name,'_',recon{1,2}{indx}]));
                                end
                                
                            end
                        elseif strcmp(pipelines{folders},'preprocessing')
                            anat_dir = fullfile(new_path,pipelines{folders},'anat');
                            log_dir = fullfile(new_path,pipelines{folders},'log');
                            if ~exist(anat_dir,'dir')
                                mkdir(anat_dir)
                                
                            end
                            if ~exist(log_dir,'dir')
                                mkdir(log_dir)
                            end
                            
                            if ismember(files_to_move{files},preproc{:,1})
                                %corresponding index of the new pat
                                indx = cellfun(@(x)strcmp(x,files_to_move{files}),preproc{:,1});
                                if endsWith(files_to_move{files},'.nii')
                                    %first move%
                                    if exist(fullfile(source_patient,files_to_move{files}),'file')
                                        if move
                                           movefile(fullfile(source_patient,files_to_move{files}),fullfile(anat_dir,files_to_move{files}));
                                        else
                                            copyfile(fullfile(source_patient,files_to_move{files}),anat_dir)
                                        end
                                        %then rename%
                                        
                                        if flag
                                            movefile(fullfile(new_path,pipelines{folders},'anat',files_to_move{files}),fullfile(anat_dir,['sub-',patient_name,'_',preproc{1,2}{indx}]));
                                        else
                                            movefile(fullfile(new_path,pipelines{folders},'anat',files_to_move{files}),fullfile(anat_dir,[patient_name,'_',preproc{1,2}{indx}]));
                                        end
                                    end
                                elseif endsWith(files_to_move{files},'.log')
                                    if exist(fullfile(source_patient,files_to_move{files}),'file')
                                        if move
                                            movefile(fullfile(source_patient,files_to_move{files}),fullfile(log_dir,files_to_move{files}));
                                        else
                                            copyfile(fullfile(source_patient,files_to_move{files}),log_dir)
                                        end
                                        if flag
                                            movefile(fullfile(new_path,pipelines{folders},'log',files_to_move{files}),fullfile(log_dir,['sub-',patient_name,'_',preproc{1,2}{indx}]));
                                        else
                                            movefile(fullfile(new_path,pipelines{folders},'log',files_to_move{files}),fullfile(log_dir,[patient_name,'_',preproc{1,2}{indx}]));
                                        end
                                    end
                                end
                            end
                        elseif strcmp(pipelines{folders},'prefs')
                            if ismember(files_to_move{files},prefs{:,1})
                          
                                %corresponding index of the new pat
                                indx = cellfun(@(x)strcmp(x,files_to_move{files}),prefs{:,1});
                                if exist(fullfile(source_patient,files_to_move{files}),'file')
                                    if move
                                        movefile(fullfile(source_patient,files_to_move{files}),fullfile(new_path,pipelines{folders},files_to_move{files}));
                                    else
                                        copyfile(fullfile(source_patient,files_to_move{files}),fullfile(new_path,pipelines{folders}));
                                    end
                                    if flag
                                        movefile(fullfile(new_path,pipelines{folders},files_to_move{files}),fullfile(new_path,pipelines{folders},['sub-',patient_name,'_',prefs{1,2}{indx}]));
                                    else
                                        movefile(fullfile(new_path,pipelines{folders},files_to_move{files}),fullfile(new_path,pipelines{folders},[patient_name,'_',prefs{1,2}{indx}]));
                                    end
                                end
                            end
                        elseif strcmp(pipelines{folders},'log')
                            if exist(fullfile(source_patient,'ea_methods.txt'),'file')
                                if move
                                    movefile(fullfile(source_patient,'ea_methods.txt'),fullfile(new_path,pipelines{folders},'ea_methods.txt'));
                                else
                                    copyfile(fullfile(source_patient,'ea_methods.txt'),fullfile(new_path,pipelines{folders}));
                                end
                                if flag
                                    movefile(fullfile(new_path,pipelines{folders},'ea_methods.txt'),fullfile(new_path,pipelines{folders},['sub-',patient_name,'_',files_to_move{files}]));
                                else
                                    movefile(fullfile(new_path,pipelines{folders},'ea_methods.txt'),fullfile(new_path,pipelines{folders},[patient_name,'_',files_to_move{files}]));
                                end
                            end
                        end
                    end
                end
            otherwise
                for j=1:length(dicom_sessions)
                    new_path = fullfile(dest_patient,subfolder_cell{subfolders},patient_name,'DICOM',dicom_sessions{j});
                    if ~exist(new_path,'dir')
                        mkdir(new_path)
                    end
                    if ~isdicom
                        disp("There are no dicom images, source data folder will be empty")
                    else
                        if move
                            movefile(dicom_source,new_path)
                        else
                            copyfile(dicom_source,new_path)
                        end
                    
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

