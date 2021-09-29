function legacy2bids(source,dest,isdicom,dicom_source,doDcmConv)
%SANTIY: first check the existence of your source and dest directory. If it is not there,
%create a new directory.
for j=1:length(source)
    if ~exist(source{j},'dir')
        warndlg("The source directory you have specified does not exist")
    else
        addpath(source{j})
    end
    if isempty(dest)
        dest{j} = source{j};
        addpath(dest{j});
    end
end
%when the app is used to run this script, there is an isempty check in the
%app also. This conflicts with whether it is empty in this script (i.e., if
%the user uses the app and doesn't input an output dir, then the output dir
%is automatically changed to the input dir. So this code will run. 
if exist('dest','var') 
   for dest_dir = 1:length(dest)
       if ~exist(dest{dest_dir},'dir')
            addpath(dest{dest_dir});
            mkdir(dest{dest_dir});
       end
   end
end
%if you have dicom files and have not provided a dicom source directory
%(this can also just be source directory,so refactor this) then throw
%an error
if isdicom && ~exist('dicom_source','var')
    warndlg("You have specified that you want dicom import, but have not specified a dicom source file. Please specific your dicom source directory!")
end

%define names of the new directorey structure
modes = {'anat','func','dwi'};
sessions = {'ses-preop','ses-postop'};
if exist('doDcmConv','var') && doDcmConv
    subfolder_cell = {'sourcedata','legacy_rawdata','derivatives'};
else
    subfolder_cell = {'sourcedata','rawdata','derivatives'};
end
pipelines = {'brainshift','coregistration','normalization','reconstruction','preprocessing','prefs','log','export','stimulations','headmodel','miscellaneous','ftracking'};
%mapping will allow quick reference of the files to move: also, certain
%modalities have specific bids naming.
legacy_modalities = {'t1.nii','t2.nii','pd.nii','ct.nii','tra.nii','cor.nii','sag.nii','fgatir.nii','fa.nii','dti.nii','dti.bval','dti.bvec','t2star.nii'};
bids_modalities = {'T1w','T2w','PDw','CT','ax','cor','sag','FGATIR','fa','dwi','dwi.bval','dwi.bvec','T2starw'};
rawdata_containers = containers.Map(legacy_modalities,bids_modalities);
[brainshift,coregistration,normalization,preprocessing,reconstruction,prefs,stimulations,headmodel,miscellaneous,ftracking] = create_bids_mapping();
%these files should be converted to .json
%data structure for excel sheet later on
derivatives_cell = {};
%dir without dots is a modified dir command that with return the
%directory contents without the dots. those evaluate to folders and
%cause issues.
%by default, use the sub- prefix. If the patient already has a sub-, this
%flag is set to zero.

%get all patient names, again this will be changed based on the dataset
%fetcher command
%perform moving for each pat
tic
log_path = fullfile(dest{1},'derivatives','leaddbs','logs');
if ~exist(log_path,'dir')
    mkdir(log_path)
end
for patients = 1:length(source)
    source_patient = source{patients};
    if patients <= length(dest)
        dest_patient = dest{patients};
    end
    if isdicom
        dicom_patient = dicom_source{patients};
    end
    [~,patient_name,~] = fileparts(source_patient);
    if ~startsWith(patient_name,'sub-')
       patient_name = ['sub-',erase(patient_name,'_')];
    end
    files_in_pat_folder = dir_without_dots(source_patient);
    file_names = {files_in_pat_folder.name};
    file_index = 1;
    for j=1:length(file_names)
        %add filenames. However, some files are inside folder (checkrreg
        %and scrf). They will be handled a bit later.
        files_to_move{file_index,1} = file_names{j};
        file_index = file_index + 1;
        %now let's deal with subfolders
    end
    %collect directory names inside the patient folder.
    dir_names = {files_in_pat_folder([files_in_pat_folder.isdir]).name};
    for j=1:length(dir_names)
        if j==1
            disp("Migrating Atlas folder...");
            if any(ismember(dir_names,'WarpDrive'))
                movefile(fullfile(source_patient,'WarpDrive'),fullfile(source_patient,'warpdrive'));
                disp("Migrating warpdrive folder...");
            end
            if any(ismember(dir_names,'stimulations'))
                %if mni dir exists
                if ~exist(fullfile(source_patient,'stimulations','MNI_ICBM_2009b_NLIN_ASYM'),'dir') || ~exist(fullfile(source_patient,'stimulations','native'),'dir')
                    mkdir(fullfile(dest_patient,'derivatives','leaddbs',patient_name,'stimulation','MNI152NLin2009bAsym'));
                    copyfile(fullfile(source_patient,'stimulations'),fullfile(dest_patient,'derivatives','leaddbs',patient_name,'stimulations','MNI152NLin2009bAsym'));
                else
                    copyfile(fullfile(source_patient,'stimulations'),fullfile(dest_patient,'derivatives','leaddbs',patient_name,'stimulations'));
                end
            end
            if any(ismember(dir_names,'current_headmodel'))
                if exist(fullfile(source_patient,'headmodel'),'dir')
                    movefile(fullfile(source_patient,'current_headmodel'),fullfile(source_patient,'headmodel'));
                else %there is only a current headmodel but no headmodel
                    copyfile(fullfile(source_patient,'current_headmodel'),fullfile(dest_patient,'derivatives','leaddbs',patient_name,'headmodel'))
                end
            end
        end
       
        if strcmp(dir_names{j},'scrf') 
            new_path = fullfile(dest_patient,'derivatives','leaddbs',patient_name);
            scrf_patient = fullfile(source_patient,'scrf');
            which_pipeline = pipelines{1};
            if ~exist(fullfile(new_path,which_pipeline),'dir')
                mkdir(fullfile(new_path,which_pipeline))
            end
            brainshift_path = fullfile(new_path,'brainshift');
            this_folder = dir_without_dots(fullfile(source_patient,dir_names{j}));
            files_in_folder = {this_folder.name};
            for file_in_folder=1:length(files_in_folder)
               which_file = files_in_folder{file_in_folder};
               disp(files_in_folder{file_in_folder});
               if ismember(files_in_folder{file_in_folder},brainshift{:,1})
                    indx = cellfun(@(x)strcmp(x,files_in_folder{file_in_folder}),brainshift{:,1});
                    bids_name = brainshift{1,2}{indx};
                    ea_generate_datasetDescription(brainshift_path,'derivatives')
                    derivatives_cell = move_derivatives2bids(scrf_patient,new_path,which_pipeline,which_file,patient_name,bids_name,derivatives_cell); 
               end
            end
                
        elseif ~strcmp(dir_names{j},'current_headmodel') %already handled
            if ~exist(fullfile(dest_patient,'derivatives','leaddbs',patient_name,dir_names{j}),'dir')
                copyfile(fullfile(source_patient,dir_names{j}),fullfile(dest_patient,'derivatives','leaddbs',patient_name,dir_names{j}));
            end
        end
    end
    %so we have created a list of the files to move, now we can create
    %the new dirrectory structure and move the correct files into the
    %correct "BIDS" directory
    
    ea_generate_datasetDescription(dest_patient,'root_folder')
    if ~isdicom
        generate_rawImagejson(files_to_move,patient_name,dest_patient,rawdata_containers);
    end
    %check for files without a DICOM folder: if there are anat_t1.nii
    %derivatives in the coreg folder, but no raw_anat_t1.nii - then there
    %is no way that there can a raw data dataset formed. We will
    %incorporate a simple check to ensure that in this case, the .json file
    %inside the raw data folder points to the coreg folder and not the
    %raw_anat_t1.nii.gz
   
    
    for subfolders = 1:length(subfolder_cell)
        switch subfolder_cell{subfolders}
            case 'sourcedata'
                new_path = fullfile(dest_patient,subfolder_cell{subfolders},patient_name);
                if ~exist(new_path,'dir')
                    mkdir(new_path)
                end
                if ~isdicom
                    disp("There are no dicom images, source data folder will be empty")
                else
                    disp("Copying DICOM folder...");
                    copyfile(dicom_patient,new_path)                    
                end 
            case 'derivatives'
                disp("Migrating Derivatives folder...");
                new_path = fullfile(dest_patient,subfolder_cell{subfolders},'leaddbs',patient_name);
                for files=1:length(files_to_move)
                    mod_split = strsplit(files_to_move{files},'_');
                    if ismember(mod_split,'postop')
                        session = 'ses-postop';
                    else
                        session = 'ses-preop';
                    end
                    mod_str = mod_split{end};
                    
                    mod_str = strsplit(mod_str,'.');
                    mod_str = [upper(mod_str{1}),'.nii'];
                     
                    if ismember(files_to_move{files},coregistration{:,1})
                        %corresponding index of the new pat
                        which_pipeline = pipelines{2};
                        which_file = files_to_move{files};
                        indx = cellfun(@(x)strcmp(x,files_to_move{files}),coregistration{:,1});
                        bids_name = coregistration{1,2}{indx};
                        if ~exist(fullfile(new_path,which_pipeline),'dir')
                            mkdir(fullfile(new_path,which_pipeline))
                        end
                        coreg_path = fullfile(new_path,'coregistration');
                        ea_generate_datasetDescription(coreg_path,'derivatives')
                        derivatives_cell = move_derivatives2bids(source_patient,new_path,which_pipeline,which_file,patient_name,bids_name,derivatives_cell);
                    
                    elseif ~isempty(regexp(files_to_move{files},'^coreg.*.log'))
                        derivatives_cell{end+1,1} = fullfile(source_patient,files_to_move{files});
                        derivatives_cell{end,2} = fullfile(new_path,pipelines{2},'log',files_to_move{files});
                        if ~exist(fullfile(new_path,pipelines{2},'log'),'dir')
                            mkdir(fullfile(new_path,pipelines{2},'log'));
                        end
                        if exist(fullfile(source_patient,files_to_move{files}),'file')
                            copyfile(fullfile(source_patient,files_to_move{files}),fullfile(new_path,pipelines{2},'log'));
                        end
                        
                        
                    elseif ismember(files_to_move{files},normalization{:,1})
                        %corresponding index of the new pat
                        which_file = files_to_move{files};
                        which_pipeline = pipelines{3};
                        indx = cellfun(@(x)strcmp(x,files_to_move{files}),normalization{:,1});
                        bids_name = normalization{1,2}{indx};
                        if ~exist(fullfile(new_path,which_pipeline),'dir')
                            mkdir(fullfile(new_path,which_pipeline))
                        end
                        normalization_path = fullfile(new_path,'normalization');
                        ea_generate_datasetDescription(normalization_path,'derivatives');
                        derivatives_cell = move_derivatives2bids(source_patient,new_path,which_pipeline,which_file,patient_name,bids_name,derivatives_cell);
                                %only for normalization
                    elseif ~isempty(regexp(files_to_move{files},'^normalize_.*.log'))
                        derivatives_cell{end+1,1} = fullfile(source_patient,files_to_move{files});
                        derivatives_cell{end,2} = fullfile(new_path,pipelines{3},'log',files_to_move{files});
                        if ~exist(fullfile(new_path,pipelines{3},'log'),'dir')
                            mkdir(fullfile(new_path,pipelines{3},'log'));
                        end
                        if exist(fullfile(source_patient,files_to_move{files}),'file')
                            copyfile(fullfile(source_patient,files_to_move{files}),fullfile(new_path,pipelines{3},'log'));
                        end
                            %special case for recon
                    elseif ismember(files_to_move{files},reconstruction{:,1})
                        recon_dir = fullfile(new_path,pipelines{4});
                        if ~exist(recon_dir,'dir')
                            mkdir(recon_dir)
                        end
                        ea_generate_datasetDescription(recon_dir,'derivatives')
                        %corresponding index of the new pat
                        indx = cellfun(@(x)strcmp(x,files_to_move{files}),reconstruction{:,1});
                        bids_name = reconstruction{1,2}{indx};
                        derivatives_cell{end+1,1} = fullfile(source_patient,which_file);
                        derivatives_cell{end,2} = fullfile(recon_dir,[patient_name,'_',bids_name]);
                        copyfile(fullfile(source_patient,files_to_move{files}),recon_dir)
                        movefile(fullfile(new_path,pipelines{4},files_to_move{files}),fullfile(recon_dir,[patient_name,'_',reconstruction{1,2}{indx}]));
                    
                    elseif ismember(files_to_move{files},preprocessing{:,1})
                        which_file = files_to_move{files};
                        %corresponding index of the new pat
                        indx = cellfun(@(x)strcmp(x,files_to_move{files}),preprocessing{:,1});
                        which_pipeline = pipelines{5};
                        bids_name = preprocessing{1,2}{indx};
                        if ~exist(fullfile(new_path,which_pipeline),'dir')
                            mkdir(fullfile(new_path,which_pipeline))
                        end
                        preproc_path = fullfile(new_path,'preprocessing');
                        ea_generate_datasetDescription(preproc_path,'derivatives')
                        derivatives_cell = move_derivatives2bids(source_patient,new_path,which_pipeline,which_file,patient_name,bids_name,derivatives_cell);
                        
                    
                    elseif ismember(files_to_move{files},prefs{:,1})
                        %corresponding index of the new pat
                        bids_name = [patient_name,'_','desc-','uiprefs.mat'];
                        derivatives_cell{end+1,1} = fullfile(source_patient,files_to_move{files});
                        derivatives_cell{end,2} = fullfile(new_path,pipelines{6},bids_name);
                        which_pipeline = pipelines{6};
                        if ~exist(fullfile(new_path,which_pipeline),'dir')
                            mkdir(fullfile(new_path,which_pipeline));
                        end
                        prefs_path = fullfile(new_path,'prefs');
                        ea_generate_datasetDescription(prefs_path,'derivatives')
                        copyfile(fullfile(source_patient,'ea_ui.mat'),fullfile(new_path,pipelines{6}));
                        movefile(fullfile(new_path,pipelines{6},'ea_ui.mat'),fullfile(new_path,pipelines{6},bids_name));
                            
                   
                    elseif strcmp(files_to_move{files},'ea_methods.txt') && exist(fullfile(source_patient,'ea_methods.txt'),'file')
                        bids_name = [patient_name,'_','desc-','ea_methods.txt'];
                        derivatives_cell{end+1,1} = fullfile(source_patient,files_to_move{files});
                        derivatives_cell{end,2} = fullfile(new_path,pipelines{7},bids_name);
                        which_pipeline = pipelines{7};
                        if ~exist(fullfile(new_path,which_pipeline),'dir')
                            mkdir(fullfile(new_path,which_pipeline));
                        end
                        log_path = fullfile(new_path,'log');
                        ea_generate_datasetDescription(log_path,'derivatives');
                        copyfile(fullfile(source_patient,'ea_methods.txt'),fullfile(new_path,pipelines{7}));
                        movefile(fullfile(new_path,pipelines{7},'ea_methods.txt'),fullfile(new_path,pipelines{7},bids_name));

                        
                    elseif ismember(files_to_move{files},miscellaneous{:,1})
                        derivatives_cell{end+1,1} = fullfile(source_patient,files_to_move{files});
                        derivatives_cell{end,2} = fullfile(new_path,pipelines{11},files_to_move{files});
                        which_pipeline = pipelines{11};
                        if ~exist(fullfile(new_path,which_pipeline),'dir')
                            mkdir(fullfile(new_path,which_pipeline));
                        end
                        misc_path = fullfile(new_path,'miscellaneous');
                        ea_generate_datasetDescription(misc_path,'derivatives')
                        copyfile(fullfile(source_patient,files_to_move{files}),fullfile(new_path,pipelines{11}));
                    
                    elseif ismember(files_to_move{files},ftracking{:,1})
                        which_file = files_to_move{files};
                        which_pipeline = pipelines{12};
                        indx = cellfun(@(x)strcmp(x,files_to_move{files}),ftracking{:,1});
                        bids_name = ftracking{1,2}{indx};
                        if ~exist(fullfile(new_path,which_pipeline),'dir')
                            mkdir(fullfile(new_path,which_pipeline))
                        end
                        ftracking_path = fullfile(new_path,'ftracking');
                        ea_generate_datasetDescription(ftracking_path,'derivatives')
                        derivatives_cell = move_derivatives2bids(source_patient,new_path,which_pipeline,which_file,patient_name,bids_name,derivatives_cell);
                    
                    elseif ~ismember(files_to_move{files},preprocessing{:,1}) && ~isempty(regexp(files_to_move{files},'raw_.*')) %support for other modalities in preproc
                        %other raw files go to pre-processing folder.
                        bids_name = [patient_name,'_','desc-preproc_',session,'_',mod_str];
                        if ~exist(fullfile(new_path,pipelines{5},'anat'),'dir')
                            mkdir(fullfile(new_path,pipelines{5},'anat'))
                        end
                        
                        copyfile(fullfile(source_patient,files_to_move{files}),fullfile(new_path,pipelines{5},'anat'));
                        movefile(fullfile(new_path,pipelines{5},'anat',files_to_move{files}),fullfile(new_path,pipelines{5},'anat',bids_name));
                    
                    
                    elseif ~ismember(files_to_move{files},coregistration{:,1}) && ~isempty(regexp(files_to_move{files},'^anat_.*')) %support for other modalities in coreg
                        bids_name = [patient_name,'_','space-anchorNative_desc-preproc_',session,'_',mod_str];
                        if ~exist(fullfile(new_path,pipelines{2},'anat'),'dir')
                            mkdir(fullfile(new_path,pipelines{2},'anat'))
                        end
                        copyfile(fullfile(source_patient,files_to_move{files}),fullfile(new_path,pipelines{2},'anat'));
                        movefile(fullfile(new_path,pipelines{2},'anat',files_to_move{files}),fullfile(new_path,pipelines{2},'anat',bids_name));
                   
                    elseif ~ismember(files_to_move{files},normalization{:,1}) && ~isempty(regexp(files_to_move{files},'^glanat_.*')) %support for other modalities in normalization                        
                        bids_name = [patient_name,'_','space-MNI152NLin2009bAsym_desc-preproc_',session,'_',mod_str];
                        if ~exist(fullfile(new_path,pipelines{3},'anat'),'dir')
                            mkdir(fullfile(new_path,pipelines{3},'anat'))
                        end
                        copyfile(fullfile(source_patient,files_to_move{files}),fullfile(new_path,pipelines{3},'anat'));
                        movefile(fullfile(new_path,pipelines{3},'anat',files_to_move{files}),fullfile(new_path,pipelines{3},'anat',bids_name));                    
                    end
                end
                for folders = 1:length(pipelines)
                    if strcmp(pipelines{folders},'stimulations')
                        stim_path = fullfile(new_path,'stimulations');
                         ea_generate_datasetDescription(stim_path,'derivatives');
                        %the stimulations folder should already be
                        %there in the dest directory.
                        if exist(fullfile(source_patient,pipelines{folders}),'dir') && exist(fullfile(new_path,pipelines{folders}),'dir')
                            pipeline = pipelines{folders};
                            try
                                [mni_files,native_files,derivatives_cell] = vta_walkpath(source_patient,new_path,pipeline,derivatives_cell);
                                move_mni2bids(mni_files,native_files,stimulations,'',pipeline,patient_name);
                            catch
                                disp("Your stimulation folder might be empty..."); 
                            end
                        end
                        
                    elseif strcmp(pipelines{folders},'headmodel')
                        headmodel_path = fullfile(new_path,'headmodel');
                        ea_generate_datasetDescription(headmodel_path,'derivatives');
                        if exist(fullfile(source_patient,pipelines{folders}),'dir') && exist(fullfile(new_path,pipelines{folders}),'dir')
                            headmodel_contents = dir_without_dots(fullfile(new_path,pipelines{folders}));
                            headmodel_files = {headmodel_contents.name};
                            for headmodel_file = 1:length(headmodel_files)
                                if ismember(headmodel_files{headmodel_file},headmodel{:,1})
                                    indx = cellfun(@(x)strcmp(x,headmodel_files{headmodel_file}),headmodel{:,1});
                                    derivatives_cell{end+1,1} = fullfile(source_patient,pipelines{folders},headmodel_files{headmodel_file});
                                    derivatives_cell{end,2} = fullfile(new_path,pipelines{folders},[patient_name,'_',headmodel{1,2}{indx}]);
                                    movefile(fullfile(new_path,pipelines{folders},headmodel_files{headmodel_file}),fullfile(new_path,pipelines{folders},[patient_name,'_',headmodel{1,2}{indx}]));
                                end
                            end
                        
                            curr_pipeline = 'current_headmodel';
                            new_path = fullfile(new_path,'headmodel');
                            try
                                [mni_files,native_files,derivatives_cell] = vta_walkpath(source_patient,new_path,curr_pipeline,derivatives_cell);
                                move_mni2bids(mni_files,native_files,'',headmodel,curr_pipeline,patient_name);
                            catch
                                disp("Could not detect headmodel files...");
                            end
                        elseif exist(fullfile(source_patient,'current_headmodel'),'dir') && exist(fullfile(new_path,pipelines{folders}),'dir')
                            headmodel_contents = dir_without_dots(fullfile(new_path,pipelines{folders},'MNI_ICBM_2009b_NLIN_ASYM'));
                            headmodel_files = {headmodel_contents.name};
                            for headmodel_file = 1:length(headmodel_files)
                                if ismember(headmodel_files{headmodel_file},headmodel{:,1})
                                    indx = cellfun(@(x)strcmp(x,headmodel_files{headmodel_file}),headmodel{:,1});
                                    derivatives_cell{end+1,1} = fullfile(source_patient,pipelines{folders},headmodel_files{headmodel_file});
                                    derivatives_cell{end,2} = fullfile(new_path,pipelines{folders},[patient_name,'_',headmodel{1,2}{indx}]);
                                    movefile(fullfile(new_path,pipelines{folders},'MNI_ICBM_2009b_NLIN_ASYM',headmodel_files{headmodel_file}),fullfile(new_path,pipelines{folders},'MNI_ICBM_2009b_NLIN_ASYM',[patient_name,'_',headmodel{1,2}{indx}]));
                                end
                            end
                            
                        end
                    end
                end
            otherwise
                for i= 1:length(modes)
                    for j=1:length(sessions)
                        new_path = fullfile(dest_patient,subfolder_cell{subfolders},patient_name,sessions{j},modes{i});
                        if ~exist(new_path,'dir')
                            mkdir(new_path)
                        end
                        if strcmp(modes{i},'anat') && strcmp(sessions{j},'ses-preop')
                            disp("Migrating pre operative session data...")
                            for files=1:length(files_to_move)
                                %files to be moved into pre-op:raw_anat_*.nii
                                if regexp(files_to_move{files},'raw_anat_.*.nii')
                                    if exist(fullfile(source_patient,files_to_move{files}),'file')
                                        modality_str = strsplit(files_to_move{files},'_');
                                        modality_str = modality_str{end};
                                        try
                                            bids_name = [sessions{j},'_',rawdata_containers(modality_str),'.nii.gz'];
                                        catch
                                           modality_str = strsplit(modality_str,'.');
                                           modality_str = upper(modality_str{1});
                                           bids_name = [sessions{j},'_',modality_str,'.nii.gz'];
                                        end
                                        source_patient_path = source_patient;
                                        which_file = files_to_move{files};
                                        derivatives_cell{end+1,1} = fullfile(source_patient_path,which_file);
                                        derivatives_cell{end,2} = fullfile(new_path,[patient_name,'_',bids_name]);
                                        move_raw2bids(source_patient_path,new_path,which_file,patient_name,bids_name)
                                    end
                                end
                            end
                        elseif strcmp(modes{i},'anat') && strcmp(sessions{j},'ses-postop')
                            disp("Migrating post operative session data...")
                            for files=1:length(files_to_move)
                                if ~isempty(regexp(files_to_move{files},'raw_postop_.*.nii')) || strcmp(files_to_move{files},'postop_ct.nii')
                                    if exist(fullfile(source_patient,files_to_move{files}),'file')
                                        modality_str = strsplit(files_to_move{files},'_');
                                        modality_str = modality_str{end};
                                        try
                                            bids_name = [sessions{j},'_',rawdata_containers(modality_str),'.nii.gz'];
                                        catch
                                           modality_str = strsplit(modality_str,'.');
                                           modality_str = upper(modality_str{1});
                                           bids_name = [sessions{j},'_',modality_str,'.nii.gz'];
                                        end
                                        source_patient_path = source_patient;
                                        which_file = files_to_move{files};
                                        derivatives_cell{end+1,1} = fullfile(source_patient_path,which_file);
                                        derivatives_cell{end,2} = fullfile(new_path,[patient_name,'_',bids_name]);
                                        move_raw2bids(source_patient_path,new_path,which_file,patient_name,bids_name)
                                    end
                                end
                            end
                        
                        elseif strcmp(modes{i},'dwi') && strcmp(sessions{j},'ses-preop')
                            disp("Migrating dwi data...")
                            for files = 1:length(files_to_move)
                                if ~isempty(regexp(files_to_move{files},'^dti.[bval,bvec,nii]'))
                                    if exist(fullfile(source_patient,files_to_move{files}),'file')
                                        modality_str = strsplit(files_to_move{files},'_');
                                        modality_str = modality_str{end};
                                        try
                                            bids_name = [sessions{j},'_',rawdata_containers(modality_str)];
                                        catch
                                           modality_str = strsplit(modality_str,'.');
                                           modality_str = upper(modality_str{1});
                                           bids_name = [sessions{j},'_',modality_str,'.nii.gz'];
                                        end
                                        source_patient_path = source_patient;
                                        which_file = files_to_move{files};
                                        derivatives_cell{end+1,1} = fullfile(source_patient_path,which_file);
                                        derivatives_cell{end,2} = fullfile(new_path,[patient_name,'_',bids_name]);
                                        move_raw2bids(source_patient_path,new_path,which_file,patient_name,bids_name)
                                       
                                    end
                                end
                            end
                        end
                        
                    end
                end
        end
    end        
    disp(['Process finished for Patient:' patient_name]);
    disp("Generating excel sheet for the conversion...");
    writecell(derivatives_cell,fullfile(dest_patient,'derivatives','leaddbs',patient_name,'legacy2bids_naming.xlsx'))
    disp(['Report saved at:' fullfile(dest_patient,'derivatives','leaddbs',patient_name,'legacy2bids_naming.xlsx')]);
end
toc;

function derivatives_cell = move_derivatives2bids(source_patient_path,new_path,which_pipeline,which_file,patient_name,bids_name,derivatives_cell)
    
    %if strcmp(which_pipeline,'brainshift')
    %    source_patient_path = fullfile(source_patient_path,'scrf');
    if strcmp(which_pipeline,'coregistration')
        brainshift_log_dir = fullfile(new_path,'brainshift','log');
    end
        
    anat_dir = fullfile(new_path,which_pipeline,'anat');
    log_dir = fullfile(new_path,which_pipeline,'log');
    
    checkreg_dir = fullfile(new_path,which_pipeline,'checkreg');    
    transformations_dir = fullfile(new_path,which_pipeline,'transformations');
    if endsWith(which_file,'.nii')
        if ~exist(anat_dir,'dir')
            mkdir(anat_dir)
        end
        old_path = fullfile(source_patient_path,which_file);
        new_path = anat_dir;
        
        if exist(old_path,'file')
            %first move%
            copyfile(old_path,new_path);
            %then rename%
            disp(['Renaming file:' which_file ' to:' bids_name]);
            rename_path = fullfile(new_path,which_file);
            derivatives_cell{end+1,1} = fullfile(old_path);
            derivatives_cell{end,2} = fullfile(new_path,[patient_name,'_',bids_name]);
            movefile(rename_path,fullfile(new_path,[patient_name,'_',bids_name]));
        end
        
    elseif endsWith(which_file,'.png')
        if ~exist(checkreg_dir,'dir')
            mkdir(checkreg_dir)
        end
        old_path = fullfile(source_patient_path,which_file);
        old_path_scrf = fullfile(source_patient_path,'checkreg',which_file);
        new_path = checkreg_dir;
       
        if exist(old_path,'file')
            %first move%
            copyfile(old_path,new_path);
            %then rename%
            disp(['Renaming file:' which_file ' to:' bids_name]);
            rename_path = fullfile(new_path,which_file);
            derivatives_cell{end+1,1} = fullfile(old_path);
            derivatives_cell{end,2} = fullfile(new_path,[patient_name,'_',bids_name]);
            movefile(rename_path,fullfile(new_path,[patient_name,'_',bids_name]));
        
        elseif exist(old_path_scrf,'file')
            copyfile(old_path_scrf,new_path);
            rename_path = fullfile(new_path,which_file);
            disp(['Renaming file:' which_file ' to:' bids_name]);
            derivatives_cell{end+1,1} = fullfile(old_path);
            derivatives_cell{end,2} = fullfile(new_path,[patient_name,'_',bids_name]);
            movefile(rename_path,fullfile(new_path,[patient_name,'_',bids_name]));
            
        end
        
    elseif endsWith(which_file,'.log') || ~isempty(regexp(which_file,'.*_approved||.*_applied.mat')) || endsWith(which_file,'.txt')
        if ~exist(log_dir,'dir')
            mkdir(log_dir)
        end
        if strcmp(fullfile(source_patient_path,which_file),fullfile(source_patient_path,'ea_methods.txt'))
            bids_name = 'desc-brainshiftmethod.txt';
        elseif strcmp(fullfile(source_patient_path,which_file),fullfile(source_patient_path,'ea_coreg_approved.mat'))
            coreg_mat = load(fullfile(source_patient_path,'ea_coreg_approved.mat'));
            if isfield(coreg_mat,'brainshift')
                brainshift_method.brainshift = coreg_mat.brainshift;
                if ~exist(brainshift_log_dir,'dir')
                    mkdir(brainshift_log_dir)
                end
                save(fullfile(brainshift_log_dir,[patient_name,'_','desc-brainshiftmethod.json']),'brainshift_method')
            end
        end        
        fname_in = fullfile(source_patient_path,which_file);
        [~,fname,ext] = fileparts(bids_name);
        fname_out = fullfile(log_dir,[patient_name,'_',fname '.json']);
        file2json(fname_in,fname_out); %under leaddbs/helpers/file2json
        
    elseif endsWith(which_file,'.mat') || endsWith(which_file,'.gz') || endsWith(which_file,'.h5')
        if ~exist(transformations_dir,'dir')
            mkdir(transformations_dir)
        end
        old_path = fullfile(source_patient_path,which_file);
        new_path = transformations_dir;
        if exist(old_path,'file')
            %first move%
            copyfile(old_path,new_path);
            %then rename%
            disp(['Renaming file:' which_file ' to:' bids_name]);
            rename_path = fullfile(new_path,which_file);
            derivatives_cell{end+1,1} = fullfile(old_path);
            derivatives_cell{end,2} = fullfile(new_path,[patient_name,'_',bids_name]);
            movefile(rename_path,fullfile(new_path,[patient_name,'_',bids_name]));
            
        end
    end
    
    return 
    
function move_raw2bids(source_patient_path,new_path,which_file,patient_name,bids_name)
    if exist(fullfile(source_patient_path,which_file),'file')
        copyfile(fullfile(source_patient_path,which_file),new_path);
         if endsWith(which_file,'.nii')
            gzip(fullfile(new_path,which_file))
            ea_delete(fullfile(new_path,which_file))
            which_file = [which_file,'.gz'];
         end
         
        movefile(fullfile(new_path,which_file),fullfile(new_path,[patient_name,'_',bids_name]));
        
    end
function move_mni2bids(mni_files,native_files,stimulations,headmodel,which_pipeline,patient_name)
    if strcmp(which_pipeline,'stimulations')
        if ~isempty(mni_files)
            for mni_file = 1:length(mni_files)
                for mni_subfile = 1:length(mni_files{1,mni_file})
                    [filepath,mni_filename,ext] = fileparts(mni_files{1,mni_file}{1,mni_subfile});
                    if ismember([mni_filename,ext],stimulations{:,1})
                        indx = cellfun(@(x)strcmp(x,[mni_filename,ext]),stimulations{:,1});
                        movefile(mni_files{1,mni_file}{1,mni_subfile},fullfile(filepath,[patient_name,'_',stimulations{1,2}{indx}]));
                    end
                end
            end
        end
        if ~isempty(native_files)
            for native_file = 1:length(native_files)
                for native_subfile = 1:length(native_files{1,native_file})
                    [filepath,native_filename,ext] = fileparts(native_files{1,native_file}{1,native_subfile});
                    if ismember([native_filename,ext],stimulations{:,1})
                        indx = cellfun(@(x)strcmp(x,[native_filename,ext]),stimulations{:,1});
                        movefile(native_files{1,native_file}{1,native_subfile},fullfile(filepath,[patient_name,'_',stimulations{1,2}{indx}]));
                    end
                end
            end
        end
    elseif strcmp(which_pipeline,'current_headmodel')
        if ~isempty(mni_files)
            for mni_file = 1:length(mni_files)
                [filepath,mni_filename,ext] = fileparts(mni_files{mni_file});
                if ismember([mni_filename,ext],headmodel{:,1})
                    indx = cellfun(@(x)strcmp(x,[mni_filename,ext]),headmodel{:,1});
                    movefile(mni_files{mni_file},fullfile(filepath,[patient_name,'_',headmodel{1,2}{indx}]));
                end
            end
        end
        if ~isempty(native_files)
            for native_file = 1:length(native_files)
                [filepath,native_filename,ext] = fileparts(native_files{native_file});
                if ismember([native_filename,ext],headmodel{:,1})
                    indx = cellfun(@(x)strcmp(x,[native_filename,ext]),headmodel{:,1});
                    movefile(native_files{native_file},fullfile(filepath,[patient_name,'_',headmodel{1,2}{indx}]));
                end
            end
        end
    end



    
function generate_rawImagejson(files_to_move,patient_name,dest_patient,rawdata_containers)
    output_dir = fullfile(dest_patient,'derivatives','leaddbs',patient_name,'prefs');
    if ~exist(output_dir,'dir')
        mkdir(output_dir)
    end
    filename_out = fullfile(dest_patient,'derivatives','leaddbs',patient_name,'prefs',[patient_name,'_','desc-rawimages.json']);
    fout_fid = fopen(filename_out,'w');
    %special_case
    if ~isempty(regexp(files_to_move,'anat_t1.nii')) && all(cellfun('isempty',regexp(files_to_move,'raw_anat_t1.nii')))
        anat_files_selected.preop.T1w = [patient_name,'_','space-anchorNative_desc-preproc_ses-preop_T1w.nii'];
    end
    %collect all the raw files from the files to move.
    for i=1:length(files_to_move)
        if ~isempty(regexp(files_to_move{i},'raw_anat_.*.nii')) 
            sessions = 'ses-preop';
            modality_str = strsplit(files_to_move{i},'_');
            modality_str = modality_str{end};
            try
                bids_name = [patient_name,'_',sessions,'_',rawdata_containers(modality_str)];
                rawdata_fieldname = strsplit(rawdata_containers(modality_str),'.');
                rawdata_fieldname = rawdata_fieldname{1};
            catch
                modality_str = strsplit(modality_str,'.');
                modality_str = upper(modality_str{1});
                bids_name = [patient_name,'_',sessions,'_',modality_str];
                rawdata_fieldname = modality_str;
            end
            
            anat_files_selected.preop.anat.(rawdata_fieldname) = bids_name;
        elseif ~isempty(regexp(files_to_move{i},'raw_postop_.*')) || strcmp(files_to_move{i},'postop_ct.nii')
            sessions = 'ses-postop';
            modality_str = strsplit(files_to_move{i},'_');
            modality_str = modality_str{end};
            try
                bids_name = [patient_name,'_',sessions,'_',rawdata_containers(modality_str)];
                rawdata_fieldname = strsplit(rawdata_containers(modality_str),'.');
                rawdata_fieldname = rawdata_fieldname{1};
            catch
                modality_str = strsplit(modality_str,'.');
                modality_str = upper(modality_str);
                bids_name = [patient_name,'_',sessions,'_',modality_str];
                rawdata_fieldname = modality_str{1};
            end            
            anat_files_selected.postop.anat.(rawdata_fieldname) = bids_name;
        end
        
    end
    encodejson = jsonencode(anat_files_selected);
    fprintf(fout_fid,encodejson);
        
    
 


        



