function ea_legacy2bids(source,dest,isdicom,dicom_source,doDcmConv,doOnlyRaw)


%SANTIY: first check the existence of your source and dest directory. If it is not there,
%create a new directory.
for j=1:length(source)
    if ~exist(source{j},'dir')
        warndlg("The source directory you have specified does not exist")
    else
        addpath(source{j})
    end
end


%SANITY for dest directories
if ~exist('dest','var')
    disp("Please select output directory")
    return
else
   if iscell(dest)
       dest = dest{1};
   end
   if ~isfolder(dest)
       addpath(dest);
       mkdir(dest);
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
if  doDcmConv && ~doOnlyRaw
    subfolder_cell = {'sourcedata','legacy_rawdata','derivatives'};
elseif doOnlyRaw && ~doDcmConv
    subfolder_cell = {'rawdata'};
elseif  doOnlyRaw && doDcmConv
    subfolder_cell = {'legacy_rawdata'};
elseif ~doOnlyRaw && ~doDcmConv
    subfolder_cell = {'rawdata','derivatives'};
end
pipelines = {'brainshift','coregistration','normalization','reconstruction','preprocessing','prefs','log','export','stimulations','headmodel','miscellaneous','ftracking'};


%mapping will allow quick reference of the files to move: also, certain
%modalities have specific bids naming.
legacy_modalities = {'t1.nii','t2.nii','pd.nii','ct.nii','tra.nii','cor.nii','sag.nii','fgatir.nii','fa.nii','dti.nii','dti.bval','dti.bvec','t2star.nii'};
bids_modalities = {'T1w','T2w','PDw','CT','acq-ax_MRI','acq-cor_MRI','acq-sag_MRI','FGATIR','fa','dwi','dwi.bval','dwi.bvec','T2starw'};
rawdata_containers = containers.Map(legacy_modalities,bids_modalities);
[brainshift,coregistration,normalization,preprocessing,reconstruction,prefs,stimulations,headmodel,miscellaneous,ftracking] = ea_create_bids_mapping();

%data structure for excel sheet later on
derivatives_cell = {};

%dir without dots is a modified dir command that with return the
%directory contents without the dots. those evaluate to folders and
%cause issues.
%by default, use the sub- prefix. If the patient already has a sub-, this
%flag is set to zero.


%perform moving for each pat
tic
log_path = fullfile(dest,'derivatives','leaddbs','logs');
if ~exist(log_path,'dir')
    mkdir(log_path)
end


for patients = 1:length(source)
    %source patient filepath
    source_patient = source{patients};
    
    %dest directory is already specified
    if isdicom
        dicom_patient = dicom_source{patients};
    end
    
    %get and refactor patient names. specially, don't allow
    %'_' or '-'
    [~,patient_name,~] = fileparts(source_patient);
    if ~startsWith(patient_name,'sub-')
        patient_name = ['sub-', regexprep(patient_name, '[\W_]', '')];
    end
    spaces_in_pat_name = isspace(patient_name);
    patient_name = patient_name(~spaces_in_pat_name);
    disp(['Processing patient: ' patient_name]);
    %handle the files in the patient folder (directories are handled later)
    %creates a cell of the files to move, later, we can create
    %the new dirrectory structure and move the files into the
    %correct "BIDS" directory

    files_to_move = {};
    files_in_pat_folder = dir_without_dots(source_patient);
    file_names = {files_in_pat_folder.name};
    file_index = 1;
    for j=1:length(file_names)
        if ~isdir(fullfile(source_patient,file_names{j}))
            files_to_move{file_index,1} = file_names{j};
            file_index = file_index + 1;
        end
    end
    
    %collect directories inside the patient folder.
    dir_names = {files_in_pat_folder([files_in_pat_folder.isdir]).name};
    for j=1:length(dir_names)
        if j==1
            
            if any(ismember(dir_names,'WarpDrive'))
                disp("Migrating warpdrive folder...");
                movefile(fullfile(source_patient,'WarpDrive'),fullfile(dest,'derivatives','leaddbs',patient_name,'warpdrive'))
            end
            if any(ismember(dir_names,'stimulations'))
                %if mni dir exists
                if ~exist(fullfile(source_patient,'stimulations','MNI_ICBM_2009b_NLIN_ASYM'),'dir') || ~exist(fullfile(source_patient,'stimulations','native'),'dir')
                    mkdir(fullfile(dest,'derivatives','leaddbs',patient_name,'stimulations','MNI152NLin2009bAsym'));
                    copyfile(fullfile(source_patient,'stimulations'),fullfile(dest,'derivatives','leaddbs',patient_name,'stimulations','MNI152NLin2009bAsym'));
                else
                    copyfile(fullfile(source_patient,'stimulations'),fullfile(dest,'derivatives','leaddbs',patient_name,'stimulations'));
                end
            end
            if any(ismember(dir_names,'current_headmodel'))
                if exist(fullfile(source_patient,'headmodel'),'dir')
                    movefile(fullfile(source_patient,'current_headmodel'),fullfile(source_patient,'headmodel'));
                else %there is only a current headmodel but no headmodel
                    copyfile(fullfile(source_patient,'current_headmodel'),fullfile(dest,'derivatives','leaddbs',patient_name,'headmodel'))
                end
            end
            if any(ismember(dir_names,'DICOM'))
              if ~exist(fullfile(dest,'sourcedata',patient_name,'DICOM'),'dir')
                  mkdir(fullfile(dest,'sourcedata',patient_name,'DICOM'))
              end
              copyfile(fullfile(source_patient,'DICOM'),fullfile(dest,'sourcedata',patient_name,'DICOM'))
            end
        end
        %handle brainshift copy and rename
        if strcmp(dir_names{j},'scrf') 
            new_path = fullfile(dest,'derivatives','leaddbs',patient_name);
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
               if ismember(files_in_folder{file_in_folder},brainshift{:,1})
                    indx = cellfun(@(x)strcmp(x,files_in_folder{file_in_folder}),brainshift{:,1});
                    bids_name = brainshift{1,2}{indx};
                    derivatives_cell = move_derivatives2bids(scrf_patient,new_path,which_pipeline,which_file,patient_name,bids_name,derivatives_cell); 
               end
            end
        %other directories the user may have    
        else
            if ~strcmp(dir_names{j},'current_headmodel') && ~strcmp(dir_names{j},'DICOM') && ~strcmp(dir_names{j},'checkreg')%already handled
                if ~exist(fullfile(dest,'derivatives','leaddbs',patient_name,dir_names{j}),'dir')
                    copyfile(fullfile(source_patient,dir_names{j}),fullfile(dest,'derivatives','leaddbs',patient_name,dir_names{j}));
                end
            end
        end
    end
    
    %generate the dataset description in the root_folder
    ea_generate_datasetDescription(dest,'root_folder')
    if ~isdicom && ~doOnlyRaw
        %generate raw image json in the raw data folder
        generate_rawImagejson(files_to_move,patient_name,dest,rawdata_containers,doOnlyRaw);
    elseif ~isdicom && doOnlyRaw
        generate_rawImagejson(files_to_move,patient_name,dest,rawdata_containers,doOnlyRaw);
    end
 
    for subfolders = 1:length(subfolder_cell)
        switch subfolder_cell{subfolders}
            case 'sourcedata'
                new_path = fullfile(dest,subfolder_cell{subfolders},patient_name);
                if ~exist(new_path,'dir')
                    mkdir(new_path)
                end
                if ~isdicom
                    disp("There are no dicom images, source data folder will be empty")
                else
                    if exist(dicom_patient, 'dir') 
                        disp("Copying DICOM folder...");
                        copyfile(dicom_patient,new_path)
                    elseif exist(fullfile(source_patient,'dicom'),'dir')
                        disp("Copying dicom folder...");
                        copyfile(fullfile(source_patient,'dicom'),new_path)
                    end
                end 
            case 'derivatives'
                disp("Migrating Derivatives folder...");
                new_path = fullfile(dest,subfolder_cell{subfolders},'leaddbs',patient_name);
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
                    
                    %coregistration
                    if ismember(files_to_move{files},coregistration{:,1})
                        
                        which_pipeline = pipelines{2};
                        which_file = files_to_move{files};
                        indx = cellfun(@(x)strcmp(x,files_to_move{files}),coregistration{:,1});
                        bids_name = coregistration{1,2}{indx};
                        if ~exist(fullfile(new_path,which_pipeline),'dir')
                            mkdir(fullfile(new_path,which_pipeline))
                        end
                        
                        derivatives_cell = move_derivatives2bids(source_patient,new_path,which_pipeline,which_file,patient_name,bids_name,derivatives_cell);
                    %coregistration: log, no fixed naming pattern and hence
                    %in an elseif command
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
                        
                        copyfile(fullfile(source_patient,'ea_methods.txt'),fullfile(new_path,pipelines{7}));
                        movefile(fullfile(new_path,pipelines{7},'ea_methods.txt'),fullfile(new_path,pipelines{7},bids_name));

                        
                    elseif ismember(files_to_move{files},miscellaneous{:,1})
                        derivatives_cell{end+1,1} = fullfile(source_patient,files_to_move{files});
                        derivatives_cell{end,2} = fullfile(new_path,pipelines{11},files_to_move{files});
                        which_pipeline = pipelines{11};
                        if ~exist(fullfile(new_path,which_pipeline),'dir')
                            mkdir(fullfile(new_path,which_pipeline));
                        end
                        
                        copyfile(fullfile(source_patient,files_to_move{files}),fullfile(new_path,pipelines{11}));
                    
                    elseif ismember(files_to_move{files},ftracking{:,1})
                        which_file = files_to_move{files};
                        which_pipeline = pipelines{12};
                        indx = cellfun(@(x)strcmp(x,files_to_move{files}),ftracking{:,1});
                        bids_name = ftracking{1,2}{indx};
                        if ~exist(fullfile(new_path,which_pipeline),'dir')
                            mkdir(fullfile(new_path,which_pipeline))
                        end
                        
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
                       
                        %the stimulations folder should already be
                        %there in the dest directory.
                        if exist(fullfile(source_patient,pipelines{folders}),'dir') && exist(fullfile(new_path,pipelines{folders}),'dir')
                            pipeline = pipelines{folders};
                            try
                                [mni_files,native_files,derivatives_cell] = ea_vta_walkpath(source_patient,new_path,pipeline,derivatives_cell);
                                move_mni2bids(mni_files,native_files,stimulations,'',pipeline,patient_name);
                            catch
                                disp("Your stimulation folder might be empty..."); 
                            end
                        end
                        
                    elseif strcmp(pipelines{folders},'headmodel')
                        
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
                                [mni_files,native_files,derivatives_cell] = ea_vta_walkpath(source_patient,new_path,curr_pipeline,derivatives_cell);
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
                if doOnlyRaw
                    preop_raw_str = '(postop||tra||sag||cor||ct)';
                    tmp_preop = cellfun('isempty', regexpi(files_to_move, preop_raw_str));
                    postop_raw_str = '(postop||tra||sag||cor||ct)';
                    tmp_postop = ~cellfun('isempty', regexpi(files_to_move, postop_raw_str));
                else
                    preop_raw_str = 'raw_anat_.*.nii';
                    tmp_preop = ~cellfun('isempty', regexpi(files_to_move, preop_raw_str));
                    postop_raw_str = '^(raw_postop|postop_).*.nii';
                    tmp_postop = ~cellfun('isempty', regexpi(files_to_move, postop_raw_str));
                end
                matching_files_preop = files_to_move(tmp_preop); %remove postop files and get only preop
                matching_files_postop = files_to_move(tmp_postop);
                for i= 1:length(modes)
                    for j=1:length(sessions)
                        new_path = fullfile(dest,subfolder_cell{subfolders},patient_name,sessions{j},modes{i});
                        if ~exist(new_path,'dir')
                            mkdir(new_path)
                        end
                        if strcmp(modes{i},'anat') && strcmp(sessions{j},'ses-preop')
                            disp("Migrating pre operative session data...")
                            %files to be moved into pre-op:raw_anat_*.nii
                            for matching_files = 1:length(matching_files_preop)
                                if exist(fullfile(source_patient,matching_files_preop{matching_files}),'file')
                                    modality_str = strsplit(matching_files_preop{matching_files},'_');
                                    modality_str = modality_str{end};
                                    try
                                        bids_name = [sessions{j},'_',rawdata_containers(modality_str),'.nii.gz'];
                                    catch
                                        modality_str = strsplit(modality_str,'.');
                                        modality_str = upper(modality_str{1});
                                        bids_name = [sessions{j},'_',modality_str,'.nii.gz'];
                                    end
                                    source_patient_path = source_patient;
                                    which_file = matching_files_preop{matching_files};
                                    derivatives_cell{end+1,1} = fullfile(source_patient_path,which_file);
                                    derivatives_cell{end,2} = fullfile(new_path,[patient_name,'_',bids_name]);
                                    move_raw2bids(source_patient_path,new_path,which_file,patient_name,bids_name)
                                end
                            end
                        elseif strcmp(modes{i},'anat') && strcmp(sessions{j},'ses-postop')
                            disp("Migrating post operative session data...")
                            for matching_files = 1:length(matching_files_postop)
                                if strcmp(matching_files_postop{matching_files},'postop_ct.nii')
                                    postop_modality = 'CT';
                                else
                                    postop_modality = 'MRI';
                                end
                                if exist(fullfile(source_patient,matching_files_postop{matching_files}),'file')
                                    modality_str = strsplit(matching_files_postop{matching_files},'_');
                                    modality_str = modality_str{end};
                                    try
                                        bids_name = [sessions{j},'_',rawdata_containers(modality_str),'.nii.gz'];
                                    catch
                                        modality_str = strsplit(modality_str,'.');
                                        modality_str = upper(modality_str{1});
                                        bids_name = [sessions{j},'_',modality_str,'.nii.gz'];
                                    end
                                    source_patient_path = source_patient;
                                    which_file = matching_files_postop{matching_files};
                                    derivatives_cell{end+1,1} = fullfile(source_patient_path,which_file);
                                    derivatives_cell{end,2} = fullfile(new_path,[patient_name,'_',bids_name]);
                                    move_raw2bids(source_patient_path,new_path,which_file,patient_name,bids_name)
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
                raw_path = fullfile(dest,subfolder_cell{subfolders});
                ea_generate_datasetDescription(raw_path,'raw',postop_modality);
        end
    end        
    disp(['Process finished for Patient:' patient_name]);
    disp("Generating excel sheet for the conversion...");
    writecell(derivatives_cell,fullfile(dest,'derivatives','leaddbs','logs','legacy2bids_naming.xlsx'))
    disp(['Report saved at:' fullfile(dest,'derivatives','leaddbs','logs','legacy2bids_naming.xlsx')]);
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
        ea_file2json(fname_in,fname_out); %under leaddbs/helpers/file2json
        
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



    
function generate_rawImagejson(files_to_move,patient_name,dest,rawdata_containers,doOnlyRaw)
    output_dir = fullfile(dest,'derivatives','leaddbs',patient_name,'prefs');
    if ~exist(output_dir,'dir')
        mkdir(output_dir)
    end
    filename_out = fullfile(dest,'derivatives','leaddbs',patient_name,'prefs',[patient_name,'_','desc-rawimages.json']);
    fout_fid = fopen(filename_out,'w');
    %special_case
    if ~isempty(regexp(files_to_move,'anat_t1.nii')) && all(cellfun('isempty',regexp(files_to_move,'raw_anat_t1.nii'))) && ~doOnlyRaw
        anat_files_selected.preop.T1w = [patient_name,'_','space-anchorNative_desc-preproc_ses-preop_T1w.nii'];
    end
    if doOnlyRaw
        preop_raw_str = '(postop||tra||sag||cor||ct)';
        tmp_preop = cellfun('isempty', regexpi(files_to_move, preop_raw_str));
        postop_raw_str = '(postop||tra||sag||cor||ct)';
        tmp_postop = ~cellfun('isempty', regexpi(files_to_move, postop_raw_str));
    else
        preop_raw_str = 'raw_anat_.*.nii';
        tmp_preop = ~cellfun('isempty', regexpi(files_to_move, preop_raw_str));
        postop_raw_str = '^(raw_postop|postop_).*.nii';
        tmp_postop = ~cellfun('isempty', regexpi(files_to_move, postop_raw_str));
    end
    matching_files_preop = files_to_move(tmp_preop); %remove postop files and get only preop
    matching_files_postop = files_to_move(tmp_postop);
    %collect all the raw files from the files to move.
    for matching_files=1:length(matching_files_preop)
        
        sessions = 'ses-preop';
        modality_str = strsplit(matching_files_preop{matching_files},'_');
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
    end
    for matching_files = 1:length(matching_files_postop)
        sessions = 'ses-postop';
        modality_str = strsplit(matching_files_postop{matching_files},'_');
        modality_str = modality_str{end};
        try
            bids_name = [patient_name,'_',sessions,'_',rawdata_containers(modality_str)];
            rawdata_fieldname = strsplit(rawdata_containers(modality_str),'.');
            rawdata_fieldname = rawdata_fieldname{1};
            if strcmp(rawdata_fieldname,'acq-ax_MRI') || strcmp(rawdata_fieldname,'acq-cor_MRI') || strcmp(rawdata_fieldname,'acq-sag_MRI')
                rawdata_fieldname = strsplit(rawdata_fieldname,'-');
                rawdata_fieldname = rawdata_fieldname{end};
                %matching_fieldname = regexp(rawdata_fieldname,'\w*ax|cor|sag\w*','match');
                %rawdata_fieldname = matching_fieldname{1};
            end
        catch
            modality_str = strsplit(modality_str,'.');
            modality_str = upper(modality_str);
            bids_name = [patient_name,'_',sessions,'_',modality_str];
            rawdata_fieldname = modality_str{1};
            if strcmp(rawdata_fieldname,'acq-ax_MRI') || strcmp(rawdata_fieldname,'acq-cor_MRI') || strcmp(rawdata_fieldname,'acq-sag_MRI')
                rawdata_fieldname = strsplit(rawdata_fieldname,'-');
                rawdata_fieldname = rawdata_fieldname{end};
                %matching_fieldname = regexp(rawdata_fieldname,'\w*ax|cor|sag\w*','match');
                %rawdata_fieldname = matching_fieldname{1};
            end
        end
        anat_files_selected.postop.anat.(rawdata_fieldname) = bids_name;
    end
    
    encodejson = jsonencode(anat_files_selected);
    fprintf(fout_fid,encodejson);
        
    
 


        



