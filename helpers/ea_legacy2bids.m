function subjId = ea_legacy2bids(source,dest,doOnlyRaw)
%%This function migrates a classic LEAD-DBS dataset, whether fully
%%processed or raw into a BIDS-STYLE dataset.The BIDSified version is
%%integral for the future releases of BIDS.
%% Parameters:
%% i)  source: full path of the source dataset (classic lead-dbs), as a cell structure. Multiple entries may be provided
%% ii) dest:   full path of the destination dataset (BIDSified lead-dbs), as a cell structure. One overarching directory should be specified
%% iii) doOnlyRaw: Boolean, 1 if your dataset only contains RAW nifti files, else 0
if ~iscell(source)
    source = {source};
end
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

%define names of the new directorey structure
pipelines = {'brainshift','coregistration','normalization','reconstruction','preprocessing','prefs','log','export','stimulations','headmodel','miscellaneous','ftracking'};
fmri_keywords = {'(Yeo 2011)','HCP 612','HCP 1200','PPMI 74','Depression 38'};

%mapping will allow quick reference of the files to move: also, certain
%modalities have specific bids naming.
legacy_modalities = {'t1','t2star','pd','ct','tra','cor','sag','fgatir','fa','dti','dti.bval','dti.bvec','t2','flair','inv','swi'};
%legacy_modalities = {'t1.nii','t2.nii','pd.nii','ct.nii','tra.nii','cor.nii','sag.nii','fgatir.nii','fa.nii','dti.nii','dti.bval','dti.bvec','t2star.nii'};
bids_modalities = {'T1w','T2starw','PDw','CT','acq-ax_MRI','acq-cor_MRI','acq-sag_MRI','FGATIR','fa','dwi','dwi.bval','dwi.bvec','T2w','FLAIR','INV','SWI'};
rawdata_containers = containers.Map(legacy_modalities,bids_modalities);
[~,brainshift,coregistration,normalization,preprocessing,reconstruction,prefs,stimulations,headmodel,ftracking,lead_mapper] = ea_create_bids_mapping();
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
ea_mkdir(log_path);

%support for lead group
subjId = cell(length(source), 1);
for patients = 1:length(source)
    %source patient filepath
    source_patient = source{patients};

    modes = {'anat','func','dwi'};
    sessions = {'ses-preop','ses-postop'};
    if doOnlyRaw
        subfolder_cell = {'rawdata'};
    else
        subfolder_cell = {'derivatives','rawdata'};
    end


    %get and refactor patient names. specially, don't allow
    %'_' or '-'
    [~,patient_name,~] = fileparts(source_patient);
    if ~startsWith(patient_name,'sub-') && ~startsWith(patient_name,'sub')
        patient_name = ['sub-', regexprep(patient_name, '[\W_]', '')];
    else
        patient_name = regexprep(patient_name, '[\W_]', '');
        patient_name = strrep(patient_name,'sub','sub-');
    end
    spaces_in_pat_name = isspace(patient_name);
    patient_name = patient_name(~spaces_in_pat_name);
    subjId{patients} = strrep(patient_name,'sub-','');
    disp(['Processing patient: ' patient_name]);
    %handle the files in the patient folder (directories are handled later)
    %creates a cell of the files to move, later, we can create
    %the new dirrectory structure and move the files into the
    %correct "BIDS" directory

    tag_cell = {}; %initializing cell for the tags
    mod_cell = {}; %initializing cell for the mods

    files_in_pat_folder = dir_without_dots(source_patient); %all the files which do not start with '.'
    files_to_move = {files_in_pat_folder(~[files_in_pat_folder.isdir]).name}';
    files_to_move = filterLegacyFiles(files_to_move);

    for j=1:length(files_to_move)
       if ~any(contains(files_to_move{j},'\w*(ct|tra|cor|sag)\w*')) %find a mapping between tags and modalities (for e.g., tag for T1w is ax, therefore tag = {'T1w.nii'}, mod = {'ax'})
           if any(regexpi(files_to_move{j},'raw_anat_.*\.nii$')) || any(regexpi(files_to_move{j},'^anat_.*\.nii$')) || doOnlyRaw  %we already know their tags in the case of cor,tra,sag
               to_match = files_to_move{j};
               bids_mod = add_mod(to_match,legacy_modalities,rawdata_containers);
               tag = ea_checkacq(fullfile(source_patient,files_to_move{j})); %function for modalities, use of fslHD
               tag_cell{end+1} = tag;
               mod_cell{end+1} = bids_mod;
           end
       end
    end
    if ~doOnlyRaw
        if any(ismember(files_to_move,'glpostop_tra.nii')) || any(ismember(files_to_move,'glpostop_sag.nii')) || any(ismember(files_to_move,'glpostop_cor.nii'))
            MRCT = 1;
        elseif any(ismember(files_to_move,'glpostop_ct.nii'))
            MRCT = 2;
        end
    end

    files_to_move = reorderfiles(files_to_move);

    %collect directories inside the patient folder.
    dir_names = {files_in_pat_folder([files_in_pat_folder.isdir]).name}; %deal with dir names
    new_path = fullfile(dest,'derivatives','leaddbs',patient_name);
    for j=1:length(dir_names)
        if strcmp(dir_names{j},'WarpDrive')
            if ~exist(fullfile(new_path,'warpdrive'),'dir')
                disp('Migrating warpdrive folder...');
                copyfile(fullfile(source_patient,'WarpDrive'),fullfile(dest,'derivatives','leaddbs',patient_name,'warpdrive'));
            end
        elseif strcmp(dir_names{j},'stimulations')
            %if mni dir exist
            if ~exist(fullfile(new_path,'stimulations'),'dir')
                copyfile(fullfile(source_patient,'stimulations'),fullfile(dest,'derivatives','leaddbs',patient_name,'stimulations'));
                if exist(fullfile(source_patient,'stimulations','MNI_ICBM_2009b_NLIN_ASYM'),'dir')
                    movefile(fullfile(dest,'derivatives','leaddbs',patient_name,'stimulations','MNI_ICBM_2009b_NLIN_ASYM'),fullfile(dest,'derivatives','leaddbs',patient_name,'stimulations','MNI152NLin2009bAsym'))
                end
            end
        elseif strcmp(dir_names{j},'current_headmodel')
                %if you have both headmodel and current_headmodel, choose
                %current_headmodel

                if exist(fullfile(source_patient,'headmodel'),'dir') %give precedence to curr_headmodel & stash headmodel into misc folder
                    misc_dir = fullfile(new_path,'miscellaneous');
                    ea_mkdir(misc_dir);
                    copyfile(fullfile(source_patient,'headmodel'),fullfile(dest,'derivatives','leaddbs',patient_name,'miscellaneous','headmodel'));
                end

                if ~exist(fullfile(new_path,'headmodel'),'dir')
                    copyfile(fullfile(source_patient,dir_names{j}),fullfile(dest,'derivatives','leaddbs',patient_name,'headmodel'));
                    if exist(fullfile(source_patient,dir_names{j},'MNI_ICBM_2009b_NLIN_ASYM'),'dir')
                        movefile(fullfile(dest,'derivatives','leaddbs',patient_name,'headmodel','MNI_ICBM_2009b_NLIN_ASYM'),fullfile(dest,'derivatives','leaddbs',patient_name,'headmodel','MNI152NLin2009bAsym'))
                    end
                end

        elseif strcmp(dir_names{j},'headmodel')
            if ~exist(fullfile(source_patient,'current_headmodel'),'dir') %if current_headmodel does not exist then this is the headmodel file
                copyfile(fullfile(source_patient,dir_names{j}),fullfile(dest,'derivatives','leaddbs',patient_name,'headmodel'));
                if exist(fullfile(source_patient,dir_names{j},'MNI_ICBM_2009b_NLIN_ASYM'),'dir')
                    movefile(fullfile(dest,'derivatives','leaddbs',patient_name,'headmodel','MNI_ICBM_2009b_NLIN_ASYM'),fullfile(dest,'derivatives','leaddbs',patient_name,'headmodel','MNI152NLin2009bAsym'))
                end
            end
        elseif strcmp(dir_names{j},'DICOM')
            ea_mkdir(fullfile(dest,'sourcedata',patient_name))
            copyfile(fullfile(source_patient,'DICOM'),fullfile(dest,'sourcedata',patient_name));
            %handle brainshift copy and rename: we already do this because
            %some of the filenames are similar to the coreg filenames and
            %in order to ensure there are no conflicts, we move scrf files
            %first.
        elseif strcmp(dir_names{j},'scrf')
            which_pipeline = pipelines{1};
            if ~exist(fullfile(new_path,which_pipeline),'dir')
                scrf_patient = fullfile(source_patient,'scrf');
                ea_mkdir(fullfile(new_path,which_pipeline));
                this_folder = dir_without_dots(fullfile(source_patient,dir_names{j}));
                files_in_folder = {this_folder.name};
                files_in_folder(contains(files_in_folder, 'rpostop_ct')) = [];
                files_in_folder(contains(files_in_folder, 'ea_methods')) = [];
                files_in_folder(ismember(files_in_folder, 'scrf.mat')) = [];
                for file_in_folder=1:length(files_in_folder)
                    which_file = files_in_folder{file_in_folder};
                    if ismember(files_in_folder{file_in_folder},brainshift{:,1})
                        indx = cellfun(@(x)strcmp(x,files_in_folder{file_in_folder}),brainshift{:,1});
                        bids_name = brainshift{1,2}{indx};
                        if MRCT == 1
                            bids_name = strrep(bids_name,'CT','MRI');
                        end
                        if contains(bids_name,'acqTag')
                            bids_name = add_tag(bids_name,mod_cell,tag_cell);
                        end
                        derivatives_cell = move_derivatives2bids(scrf_patient,new_path,which_pipeline,which_file,patient_name,bids_name,derivatives_cell);
                    else
                        misc_dir = fullfile(new_path,'miscellaneous');
                        ea_mkdir(misc_dir);
                        copyfile(fullfile(scrf_patient,files_in_folder{file_in_folder}),misc_dir)
                    end
                end
            end% delete entry from the dir names structure so as to not handle it again
            %other directories the user may have
        elseif strcmp(dir_names{j},'fiberfiltering') || strcmp(dir_names{j},'networkmapping') || strcmp(dir_names{j},'sweetspotmapping')
            ea_mkdir(fullfile(dest,'derivatives','leadgroup',patient_name,dir_names{j}));
            copyfile(fullfile(source_patient,dir_names{j}),fullfile(dest,'derivatives','leadgroup',patient_name,dir_names{j}));
          %add the .png files to the main cell which contains the files to move
        elseif strcmp(dir_names{j},'checkreg')
            checkreg_folder = dir_without_dots(fullfile(source_patient,dir_names{j}));
            files_in_checkreg_folder = {checkreg_folder.name};
            for checkreg=1:length(files_in_checkreg_folder)
                if exist('files_to_move','var')
                    files_to_move{end+1} = files_in_checkreg_folder{checkreg};
                else
                    files_to_move{checkreg} = files_in_checkreg_folder{checkreg};
                end
            end
        elseif strcmp(dir_names{j},'atlases')
            if ~exist(fullfile(dest,'derivatives','leaddbs',patient_name,dir_names{j}),'dir')
                copyfile(fullfile(source_patient,dir_names{j}),fullfile(dest,'derivatives','leaddbs',patient_name,dir_names{j}));
            end
        else
            misc_dir = fullfile(new_path,'miscellaneous');
            ea_mkdir(misc_dir);
            copyfile(fullfile(source_patient,dir_names{j}),fullfile(dest,'derivatives','leaddbs',patient_name,'miscellaneous',dir_names{j}));
        end
    end

    %generate the dataset description in the root_folder
    ea_generate_datasetDescription(dest,'root_folder')

    %handle raw,
    for subfolders = 1:length(subfolder_cell)
        switch subfolder_cell{subfolders}

            case 'derivatives'
                disp('Migrating Derivatives folder...');
                new_path = fullfile(dest,subfolder_cell{subfolders},'leaddbs',patient_name);
                for files=1:length(files_to_move)
                    which_file = files_to_move{files};
                    if contains(files_to_move{files},'postop')
                        sess_tag = 'ses-postop';
                    else
                        sess_tag = 'ses-preop';
                    end
                    %should have all the directories, even if empty to
                    %prevent errors of old datasets

                    %coregistration
                    if ismember(which_file, coregistration{:,1})
                        which_pipeline = pipelines{2};
                        indx = cellfun(@(x)strcmp(x,which_file),coregistration{:,1});
                        bids_name = coregistration{1,2}{indx};
                        if contains(bids_name,'acqTag')
                            bids_name = add_tag(bids_name,mod_cell,tag_cell);
                        end
                        %find mod of the coreg and then check if you have a
                        %raw_anat_mod in the folder.

                        derivatives_cell = move_derivatives2bids(source_patient,new_path,which_pipeline,which_file,patient_name,bids_name,derivatives_cell);

                        %coregistration: log, no fixed naming pattern and hence
                        %in an elseif command
                    elseif ~isempty(regexp(which_file, '^coreg.*\.log$', 'once'))
                        derivatives_cell{end+1,1} = fullfile(source_patient,which_file);
                        derivatives_cell{end,2} = fullfile(new_path,pipelines{2},'log',[patient_name, '_desc-',erase(which_file,'_')]);
                        if exist(fullfile(source_patient,which_file),'file')
                            ea_mkdir(fullfile(new_path,pipelines{2},'log'));
                            copyfile(derivatives_cell{end,1},derivatives_cell{end,2});
                        end

                    elseif ismember(which_file,normalization{:,1})
                        %corresponding index of the new pat
                        which_pipeline = pipelines{3};
                        indx = cellfun(@(x)strcmp(x,which_file),normalization{:,1});
                        bids_name = normalization{1,2}{indx};
                        %replace tag in the bids name
                        if contains(bids_name,'acqTag')
                            bids_name = add_tag(bids_name,mod_cell,tag_cell);
                        end
                        derivatives_cell = move_derivatives2bids(source_patient,new_path,which_pipeline,which_file,patient_name,bids_name,derivatives_cell);
                        %only for normalization
                    elseif ~isempty(regexp(which_file, '^normalize_.*\.log$', 'once'))
                        derivatives_cell{end+1,1} = fullfile(source_patient,which_file);
                        derivatives_cell{end,2} = fullfile(new_path,pipelines{3},'log',[patient_name, '_desc-',erase(which_file,'alize_')]);
                        if exist(fullfile(source_patient,which_file),'file')
                            ea_mkdir(fullfile(new_path,pipelines{3},'log'));
                            copyfile(derivatives_cell{end,1},derivatives_cell{end,2});
                        end
                        %special case for recon
                    elseif ismember(which_file,reconstruction{:,1})
                        recon_dir = fullfile(new_path,pipelines{4});
                        %corresponding index of the new pat
                        indx = cellfun(@(x)strcmp(x,which_file),reconstruction{:,1});
                        bids_name = reconstruction{1,2}{indx};
                        derivatives_cell{end+1,1} = fullfile(source_patient,which_file);
                        derivatives_cell{end,2} = fullfile(recon_dir,[patient_name,'_',bids_name]);
                        ea_mkdir(recon_dir);
                        copyfile(fullfile(source_patient,which_file),recon_dir)
                        movefile(fullfile(new_path,pipelines{4},which_file),fullfile(recon_dir,[patient_name,'_',reconstruction{1,2}{indx}]));

                    elseif ismember(which_file,preprocessing{:,1})
                        which_file = which_file;
                        %corresponding index of the new pat
                        indx = cellfun(@(x)strcmp(x,which_file),preprocessing{:,1});
                        which_pipeline = pipelines{5};
                        bids_name = preprocessing{1,2}{indx};
                        if contains(bids_name,'acqTag')
                            bids_name = add_tag(bids_name,mod_cell,tag_cell);
                        end
                        derivatives_cell = move_derivatives2bids(source_patient,new_path,which_pipeline,which_file,patient_name,bids_name,derivatives_cell);


                    elseif ismember(which_file,prefs{:,1})
                        %corresponding index of the new pat
                        bids_name = [patient_name,'_','desc-','uiprefs.mat'];
                        derivatives_cell{end+1,1} = fullfile(source_patient,which_file);
                        derivatives_cell{end,2} = fullfile(new_path,pipelines{6},bids_name);
                        which_pipeline = pipelines{6};
                        ea_mkdir(fullfile(new_path,pipelines{6}));
                        copyfile(fullfile(source_patient,'ea_ui.mat'),fullfile(new_path,pipelines{6}));
                        movefile(fullfile(new_path,pipelines{6},'ea_ui.mat'),fullfile(new_path,pipelines{6},bids_name));

                    elseif strcmp(which_file,'ea_stats.mat')
                        bids_name = [patient_name,'_','desc-','stats.mat'];
                        derivatives_cell{end+1,1} = fullfile(source_patient,which_file);
                        derivatives_cell{end,2} = fullfile(new_path,bids_name);
                        copyfile(fullfile(source_patient,'ea_stats.mat'),new_path);
                        movefile(fullfile(new_path,'ea_stats.mat'),fullfile(new_path,bids_name));

                    elseif strcmp(which_file,'ea_methods.txt') && exist(fullfile(source_patient,'ea_methods.txt'),'file')
                        bids_name = [patient_name, '_desc-methods.txt'];
                        derivatives_cell{end+1,1} = fullfile(source_patient,which_file);
                        derivatives_cell{end,2} = fullfile(new_path,pipelines{7},bids_name);
                        which_pipeline = pipelines{7};
                        ea_mkdir(fullfile(new_path,pipelines{7}));
                        copyfile(fullfile(source_patient,'ea_methods.txt'),fullfile(new_path,pipelines{7}));
                        movefile(fullfile(new_path,pipelines{7},'ea_methods.txt'),fullfile(new_path,pipelines{7},bids_name));

                    elseif ismember(which_file,ftracking{:,1})
                        which_pipeline = pipelines{12};
                        indx = cellfun(@(x)strcmp(x,which_file),ftracking{:,1});
                        bids_name = ftracking{1,2}{indx};
                        if contains(bids_name,'acqTag')
                            bids_name = add_tag(bids_name,mod_cell,tag_cell);
                        end
                        derivatives_cell = move_derivatives2bids(source_patient,new_path,which_pipeline,which_file,patient_name,bids_name,derivatives_cell);

                    elseif ~ismember(which_file,preprocessing{:,1}) && ~isempty(regexp(which_file, 'raw_.*\.nii$', 'once')) %support for other modalities in preproc
                        %other raw files go to pre-processing folder.
                        if endsWith(which_file,'.nii')
                            ext = '.nii';
                            op_dir = fullfile(new_path,pipelines{5},'anat');
                            source_path = source_patient;
                        elseif endsWith(which_file,'.png')
                            ext = '.png';
                            op_dir = fullfile(new_path,pipelines{5},'checkreg');
                            source_path = fullfile(source_patient,'checkreg');
                        end
                        ea_mkdir(op_dir);
                        bids_mod = add_mod(which_file,legacy_modalities,rawdata_containers);
                        if ~isempty(bids_mod)
                            try
                                tag = ea_checkacq(fullfile(source_path,which_file));
                                bids_name = [patient_name,'_',sess_tag,'_','desc-preproc_','acq-',tag,'_',bids_mod,ext];
                            catch
                                try_bids_name = [patient_name,'_',sess_tag,'_','desc-preproc_','acqTag','_',bids_mod,ext];
                                bids_name = add_tag(try_bids_name,mod_cell,tag_cell);
                            end
                            bids_name = CheckifAlreadyExists(op_dir,bids_name);
                            copyfile(fullfile(source_path,which_file),op_dir);
                            movefile(fullfile(op_dir,which_file),fullfile(op_dir,bids_name));
                        else
                            op_dir = fullfile(new_path,'miscellaneous');
                            ea_mkdir(op_dir);
                            copyfile(fullfile(source_path,which_file),op_dir);
                        end

                    elseif ~ismember(which_file,coregistration{:,1}) && contains(which_file, "anat_" + alphanumericsPattern)
                        if endsWith(which_file,'.nii')
                            ext = '.nii';
                            op_dir = fullfile(new_path,pipelines{2},'anat');
                            source_path = source_patient;
                        elseif endsWith(which_file,'.png')
                            ext = '.png';
                            op_dir = fullfile(new_path,pipelines{2},'checkreg');
                            source_path = fullfile(source_patient,'checkreg');
                        elseif endsWith(which_file,'.mat')
                            ext = '.mat';
                            op_dir = fullfile(new_path,pipelines{2},'transformations');
                            source_path = source_patient;
                        end
                        ea_mkdir(op_dir);
                        if endsWith(which_file,'.nii') || endsWith(which_file,'.png')
                            if endsWith(which_file,'.nii') && ~startsWith(which_file, "raw_anat" | "anat")
                                bids_mod = 'misc';
                            else
                                bids_mod = add_mod(which_file,legacy_modalities,rawdata_containers);
                            end
                            if strcmp(bids_mod,'misc')
                                op_dir = fullfile(new_path,'miscellaneous');
                                ea_mkdir(op_dir);
                                copyfile(fullfile(source_path,which_file),op_dir);
                            else
                                try
                                    tag = ea_checkacq(fullfile(source_path,which_file));
                                    bids_name = [patient_name,'_',sess_tag,'_','space-anchorNative_desc-preproc_','acq-',tag,'_',bids_mod,ext];
                                catch
                                    try_bids_name = [patient_name,'_',sess_tag,'_','space-anchorNative_desc-preproc_','acqTag','_',bids_mod,ext];
                                    bids_name = add_tag(try_bids_name,mod_cell,tag_cell);
                                end
                                bids_name = CheckifAlreadyExists(op_dir,bids_name);
                                if exist(fullfile(source_path,which_file),'file')
                                  copyfile(fullfile(source_path,which_file),op_dir);
                                elseif exist(fullfile(source_patient,which_file),'file')
                                   copyfile(fullfile(source_patient,which_file),op_dir);
                                end
                                movefile(fullfile(op_dir,which_file),fullfile(op_dir,bids_name));
                                disp(['Renaming file ' which_file ' to ' bids_name])
                                derivatives_cell{end+1,1} = fullfile(source_path,which_file);
                                derivatives_cell{end,2} = fullfile(op_dir,bids_name);

                            end
                        elseif endsWith(which_file,'.mat')
                            mat_str = regexp(which_file,'[1-9]\.mat$','split','once');
                            mat_str = mat_str{1};
                            tf = any(~cellfun('isempty',strfind(coregistration{:,1},mat_str)));
                            if tf
                                indx_arr = cellfun(@(x)strfind(x,mat_str),coregistration{:,1},'UniformOutput',false);
                                indx = find(~cellfun(@isempty,indx_arr));
                                if length(indx) == 1
                                  bids_name = [patient_name,'_',coregistration{1,2}{indx}];
                                  copyfile(fullfile(source_path,which_file),op_dir);
                                  movefile(fullfile(op_dir,which_file),fullfile(op_dir,bids_name));
                                else
                                    op_dir = fullfile(new_path,'miscellaneous');
                                    ea_mkdir(op_dir);
                                    copyfile(fullfile(source_path,which_file),op_dir);
                                end
                            end
                        end


                    elseif ~ismember(which_file,normalization{:,1}) && ~isempty(regexp(which_file, '^glanat_.*(.nii|.png)$', 'once')) %support for other modalities in normalization
                        if endsWith(which_file,'.nii')
                            ext = '.nii';
                            op_dir = fullfile(new_path,pipelines{3},'anat');
                            source_path = source_patient;
                        elseif endsWith(which_file,'.png')
                            ext = '.png';
                            op_dir = fullfile(new_path,pipelines{3},'checkreg');
                            source_path = fullfile(source_patient,'checkreg');
                        end
                        ea_mkdir(op_dir);
                        bids_mod = add_mod(which_file,legacy_modalities,rawdata_containers);
                        try
                            tag = ea_checkacq(fullfile(source_path,which_file));
                            bids_name = [patient_name,'_',sess_tag,'_','space-MNI152NLin2009bAsym_desc-preproc_','acq-',tag,'_',bids_mod,ext];
                        catch
                            try_bids_name = [patient_name,'_',sess_tag,'_','space-MNI152NLin2009bAsym_desc-preproc_','acqTag','_',bids_mod,ext];
                            bids_name = add_tag(try_bids_name,mod_cell,tag_cell);
                        end
                        bids_name = CheckifAlreadyExists(op_dir,bids_name);
                        copyfile(fullfile(source_path,which_file),op_dir);
                        movefile(fullfile(op_dir,which_file),fullfile(op_dir,bids_name));
                        %support for lead group files
                    else
                        derivatives_cell{end+1,1} = fullfile(source_patient,which_file);
                        derivatives_cell{end,2} = fullfile(new_path,pipelines{11},which_file);
                        ea_mkdir(fullfile(new_path,pipelines{11}));
                        try
                            copyfile(fullfile(source_patient,which_file),fullfile(new_path,pipelines{11}));
                        catch
                            copyfile(fullfile(source_patient,'checkreg',which_file),fullfile(new_path,pipelines{11}));
                        end
                    end
                end

                % Double check the coregmethod.json file
                modalities = checkModalities(fullfile(new_path, 'coregistration', 'anat'));
                if ~isempty(ea_regexpdir(fullfile(new_path, 'coregistration', 'anat'), '_CT\.nii$', 0, 'f'))
                    modalities = [modalities; 'CT'];
                end

                if ~isempty(modalities)
                    [~, subjPrefix] = fileparts(new_path);
                    coregMethodJson = fullfile(new_path, 'coregistration', 'log', [subjPrefix, '_desc-coregmethod.json']);
                    json = struct;
                    if isfile(coregMethodJson)
                        json = loadjson(coregMethodJson); 
                    end

                    for m=1:numel(modalities)
                        % Set approval status to 1 by default
                        if ~isfield(json, 'approval') || ~isfield(json.approval, modalities{m}) || isempty(json.approval.(modalities{m}))
                            json.approval.(modalities{m}) = 1;
                        end

                        if ~isfield(json, 'method') || ~isfield(json.method, modalities{m}) || isempty(json.method.(modalities{m}))
                            if strcmp(modalities{m}, 'CT')
                                % Fallback to default method, should never reach here
                                json.method.(modalities{m}) = 'ANTs (Avants 2008)';
                                ea_cprintf('CmdWinWarnings', 'CT coregistration method fallbacks to default: ''%s''.\n', json.method.(modalities{m}));
                            else
                                % Fallback to default method, since coregmr method is not
                                % properly stored in classic version of LeadDBS
                                json.method.(modalities{m}) = 'SPM (Friston 2007)';
                                ea_cprintf('CmdWinWarnings', 'MR coregistration method fallbacks to default: ''%s''.\n', json.method.(modalities{m}));
                            end
                        end
                    end
                    ea_mkdir(fileparts(coregMethodJson));
                    savejson('', json, struct('FileName', coregMethodJson));
                end

                for folders = 1:length(pipelines)
                    if strcmp(pipelines{folders},'stimulations')

                        %the stimulations folder should already be
                        %there in the dest directory.

                        % Fix the legacy stimulation folder without space (MNI* or native) parent folder
                        legacyStims = ea_regexpdir(fullfile(new_path,pipelines{folders}), '(?-i)^(?!(MNI|native)).+', 0, 'd');
                        if ~isempty(legacyStims)
                            templateStimFolder = ea_regexpdir(fullfile(new_path,pipelines{folders}), '(?-i)^MNI.+', 0, 'd');
                            if isempty(templateStimFolder)
                                % Create template folder in case it doesn't exist
                                templateStimFolder = {fullfile(new_path,pipelines{folders}, 'MNI152NLin2009bAsym')};
                                ea_mkdir(templateStimFolder);
                            end
                            cellfun(@(x) movefile(x, templateStimFolder{1}), legacyStims);
                        end

                        if exist(fullfile(source_patient,pipelines{folders}),'dir') && exist(fullfile(new_path,pipelines{folders}),'dir')

                            %pipeline = pipelines{folders};
                            MNI152_folders = ea_regexpdir(fullfile(new_path,pipelines{folders},'MNI152NLin2009bAsym'),'.*',0,'d');
                            native_folders = ea_regexpdir(fullfile(new_path,pipelines{folders},'native'),'.*',0,'d');
                            all_gs_folders = [MNI152_folders;native_folders];
                            for i=1:length(all_gs_folders)
                                stim_files = ea_regexpdir(all_gs_folders{i},'.*',0,'f');
                                for j=1:length(stim_files)
                                    [filepath,filename,ext] = fileparts(stim_files{j});
                                    classic_stim_file = [filename,ext];
                                    if ismember(classic_stim_file,stimulations{:,1})
                                        indx = cellfun(@(x)strcmp(x,classic_stim_file),stimulations{:,1});
                                        bids_name = [patient_name,'_',stimulations{1,2}{indx}];
                                        if contains(bids_name,'modelTag')
                                            model_name = add_model(all_gs_folders{i});
                                            bids_name = strrep(bids_name,'modelTag',model_name);
                                        end
                                        if ~isempty(indx)
                                            movefile(stim_files{j},fullfile(filepath,bids_name));
                                        end
                                    end
                                end
                                leadMapper_folder = ea_regexpdir(all_gs_folders{i},'.*',0,'d');
                                for k=1:length(leadMapper_folder)
                                    is_it_aconnectome = ea_regexpdir(leadMapper_folder{k},'.*vat_seed_compound_[df]MRI.*',0,'f');
                                    if ~isempty(is_it_aconnectome)
                                        [~,connectome_filename,ext] = fileparts(leadMapper_folder{k});
                                        try
                                            if contains(connectome_filename, fmri_keywords) % Functional connectome, subset in the name
                                                connNameSplit = strsplit(connectome_filename, '_');
                                                connectome = connNameSplit{1};
                                                subset = connNameSplit{2};
                                                bids_connectome_name = ea_getConnLabel(connectome, subset);
                                            else % Structural connectome
                                                bids_connectome_name = ea_getConnLabel(connectome_filename);
                                            end
                                            model_name = add_model(all_gs_folders{i});
                                            stimulation_folder = all_gs_folders{i};
                                            connectome_folder = fullfile(leadMapper_folder{k});
                                            %create a struct of all the properties
                                            %associated with this file
                                            tag_struct.subjID = patient_name;
                                            tag_struct.modeltag = model_name;
                                            tag_struct.conntag = bids_connectome_name;
                                            generate_bidsConnectome_name(stimulation_folder,connectome_folder,lead_mapper,stimulations,tag_struct)
                                        catch ME
                                            ea_cprintf('CmdWinWarnings', 'Connectome name "%s" may not be updated. See error report for further details.\n', connectome_filename);
                                        end
                                    end
                                end
                            end

                            %end
                        end
                    elseif strcmp(pipelines{folders},'headmodel')
                        if exist(fullfile(new_path,pipelines{folders},'MNI152NLin2009bAsym'),'dir')
                            headmodel_mni_contents = dir_without_dots(fullfile(new_path,pipelines{folders},'MNI152NLin2009bAsym'));
                            headmodel_mni_files = {headmodel_mni_contents.name};
                        else
                            headmodel_mni_files = {};
                        end
                        if exist(fullfile(new_path,pipelines{folders},'native'),'dir')
                            headmodel_native_contents = dir_without_dots(fullfile(new_path,pipelines{folders},'native'));
                            headmodel_native_files = {headmodel_native_contents.name};
                        else
                            headmodel_native_files = {};
                        end
                        which_pipeline = 'headmodel';

                        move_mni2bids(headmodel_mni_files,headmodel_native_files,'',headmodel,which_pipeline,patient_name,new_path)

                    end
                end
            otherwise
                %clean up the files to move (anything outside of nii not
                %accepted)
                for i=1:length(files_to_move)
                    if isempty(regexpi(files_to_move{i},'.*\.nii(\.gz)?$','match'))
                        files_to_move{i} = [];
                    end
                end

                files_to_move = files_to_move(~cellfun(@isempty, files_to_move));
                if doOnlyRaw
                    raw_str = '\w*(postop|ct|tra|cor|sag|auto|ignore)\w*';
                    [~,matching_files_preop] = match_exact(files_to_move,raw_str);
                    [matching_files_postop,~] = match_exact(files_to_move, raw_str);
                else
                    [matching_files_preop,~] = match_exact(files_to_move,'raw_anat_.*.nii'); %remove postop files and get only preop
                    [matching_files_postop,~] = match_exact(files_to_move,'(raw_postop_|postop_ct).*.nii');
                end

                if ~isempty(matching_files_preop) || ~isempty(matching_files_postop)
                    for i= 1:length(modes)
                        for j=1:length(sessions)
                            if strcmp(subfolder_cell{subfolders},'legacy')
                                new_path = fullfile(dest,'rawdata',subfolder_cell{subfolders},patient_name,sessions{j},modes{i});
                            else
                                new_path = fullfile(dest,subfolder_cell{subfolders},patient_name,sessions{j},modes{i});
                            end
                            if strcmp(modes{i},'anat') && strcmp(sessions{j},'ses-preop')
                                tmp_path = fullfile(new_path,'tmp');
                                ea_mkdir(tmp_path);
                                disp("Migrating pre operative session data...")
                                %files to be moved into pre-op:raw_anat_*.nii
                                for matching_files = 1:length(matching_files_preop)
                                    if exist(fullfile(source_patient,matching_files_preop{matching_files}),'file')
                                        to_match = matching_files_preop{matching_files};
                                        bids_mod = add_mod(to_match,legacy_modalities,rawdata_containers);

                                        indx = cellfun(@(x)isequal(x,bids_mod),mod_cell);
                                        unique_indx = find(indx);
                                        if length(unique_indx) > 1
                                            indx = unique_indx(1);
                                        end
                                        tag = tag_cell{indx};
                                        try_bids_name = [patient_name,'_',sessions{j},'_','acq-',tag,'_',bids_mod,'.nii.gz'];
                                        %support for multiple modalities. If a
                                        %file already exists with that name (i.e., tag & mod are the same)
                                        %we rename the old file and append a
                                        %1,2,or 3 to the acq tag of the new
                                        %file.
                                        bids_name = CheckifAlreadyExists(new_path,try_bids_name);
                                        which_file = matching_files_preop{matching_files};
                                        derivatives_cell{end+1,1} = fullfile(source_patient,which_file);
                                        derivatives_cell{end,2} = fullfile(new_path,[patient_name,'_',bids_name]);
                                        move_raw2bids(source_patient,new_path,which_file,bids_name)
                                    end
                                end
                            elseif strcmp(modes{i},'anat') && strcmp(sessions{j},'ses-postop')
                                tmp_path = fullfile(new_path,'tmp');
                                ea_mkdir(tmp_path);
                                disp("Migrating post operative session data...")
                                for matching_files = 1:length(matching_files_postop)
                                    if contains(matching_files_postop{matching_files},'ct','IgnoreCase',true)
                                        postop_modality = 'CT';
                                    else
                                        postop_modality = 'MRI';
                                    end
                                    if exist(fullfile(source_patient,matching_files_postop{matching_files}),'file')
                                        to_match = matching_files_postop{matching_files};
                                        bids_mod = add_mod(to_match,legacy_modalities,rawdata_containers);
                                        bids_name = [patient_name,'_',sessions{j},'_',bids_mod,'.nii.gz'];
                                        which_file = matching_files_postop{matching_files};
                                        derivatives_cell{end+1,1} = fullfile(source_patient,which_file);
                                        derivatives_cell{end,2} = fullfile(new_path,bids_name);
                                        move_raw2bids(source_patient,new_path,which_file,bids_name)
                                    end
                                end
                            elseif strcmp(modes{i},'dwi') && strcmp(sessions{j},'ses-preop')
                                disp("Migrating dwi data...")
                                for files = 1:length(files_to_move)
                                    if ~isempty(regexp(files_to_move{files}, '^dti\.(bval|bvec|nii)$', 'once'))
                                        if exist(fullfile(source_patient,files_to_move{files}),'file')
                                            modality_str = strsplit(files_to_move{files},'_');
                                            modality_str = lower(modality_str{end});
                                            try
                                                bids_name = [patient_name,'_',sessions{j},'_',rawdata_containers(modality_str)];
                                            catch
                                                modality_str = strsplit(modality_str,'.');
                                                modality_str = upper(modality_str{1});
                                                bids_name = [patient_name,'_',sessions{j},'_',modality_str,'.nii.gz'];
                                            end
                                            which_file = files_to_move{files};
                                            derivatives_cell{end+1,1} = fullfile(source_patient,which_file);
                                            derivatives_cell{end,2} = fullfile(new_path,bids_name);
                                            move_raw2bids(source_patient,new_path,which_file,bids_name)

                                        end
                                    end
                                end
                            end
                            if exist(tmp_path,'dir')
                                ea_delete(tmp_path);
                            end
                        end

                    end
                end

                generate_rawImagejson(patient_name,dest)
                %prev input: doDcmConv
        end
    end

    disp(['Process finished for Patient: ' patient_name]);
    disp("Generating excel sheet for the conversion...");
    writecell(derivatives_cell,fullfile(dest,'derivatives','leaddbs','logs','legacy2bids_naming.xlsx'))
    disp(['Report saved at: ' fullfile(dest,'derivatives','leaddbs','logs','legacy2bids_naming.xlsx')]);
end

toc;


function derivatives_cell = move_derivatives2bids(source_patient_path,new_path,which_pipeline,which_file,patient_name,bids_name,derivatives_cell)
    anat_dir = fullfile(new_path,which_pipeline,'anat');
    log_dir = fullfile(new_path,which_pipeline,'log');

    checkreg_dir = fullfile(new_path,which_pipeline,'checkreg');
    transformations_dir = fullfile(new_path,which_pipeline,'transformations');
    if endsWith(which_file,'.nii')
        ea_mkdir(anat_dir);
        old_path = fullfile(source_patient_path,which_file);
        new_path = anat_dir;

        if exist(old_path,'file')
            %first move%
            copyfile(old_path,new_path);
            %then rename%
            disp(['Renaming file ' which_file ' to ' bids_name]);
            rename_path = fullfile(new_path,which_file);
            derivatives_cell{end+1,1} = fullfile(old_path);
            derivatives_cell{end,2} = fullfile(new_path,[patient_name,'_',bids_name]);
            movefile(rename_path,fullfile(new_path,[patient_name,'_',bids_name]));
        end

    elseif endsWith(which_file,'.png')
        ea_mkdir(checkreg_dir);
        old_path = fullfile(source_patient_path,which_file);
        old_path_scrf = fullfile(source_patient_path,'checkreg',which_file);
        new_path = checkreg_dir;

        if exist(old_path,'file')
            %first move%
            copyfile(old_path,new_path);
            %then rename%
            disp(['Renaming file ' which_file ' to ' bids_name]);
            rename_path = fullfile(new_path,which_file);
            derivatives_cell{end+1,1} = fullfile(old_path);
            derivatives_cell{end,2} = fullfile(new_path,[patient_name,'_',bids_name]);
            movefile(rename_path,fullfile(new_path,[patient_name,'_',bids_name]));

        elseif exist(old_path_scrf,'file')
            copyfile(old_path_scrf,new_path);
            rename_path = fullfile(new_path,which_file);
            disp(['Renaming file ' which_file ' to ' bids_name]);
            derivatives_cell{end+1,1} = fullfile(old_path);
            derivatives_cell{end,2} = fullfile(new_path,[patient_name,'_',bids_name]);
            movefile(rename_path,fullfile(new_path,[patient_name,'_',bids_name]));

        end

    elseif ~isempty(regexp(which_file,'.*(_approved|_applied)\.mat$', 'once'))
        ea_mkdir(log_dir);
        % Handle the special case for branshift correction since the json
        % will be saved under brainshift/log instead of coregistration/log
        if strcmp(fullfile(source_patient_path,which_file),fullfile(source_patient_path,'ea_coreg_approved.mat'))
            coreg_mat = load(fullfile(source_patient_path,'ea_coreg_approved.mat'));
            if isfield(coreg_mat,'brainshift')
                brainshift_log_dir = fullfile(new_path,'brainshift','log');
                ea_mkdir(brainshift_log_dir);
                opt.FileName = fullfile(brainshift_log_dir,[patient_name,'_','desc-brainshiftmethod.json']);
                brainshiftLog.approval = coreg_mat.brainshift;
                % Set default method to coarse mask
                brainshiftLog.method = 'Coarse mask (Schönecker 2008)';
                savejson('',brainshiftLog,opt);
            end
        end
        fname_in = fullfile(source_patient_path,which_file);
        [~,fname,~] = fileparts(bids_name);
        fname_out = fullfile(log_dir,[patient_name,'_',fname '.json']);
        file2json(fname_in,fname_out,derivatives_cell);

    elseif endsWith(which_file,'.mat') || endsWith(which_file,'.gz') || endsWith(which_file,'.h5')
        ea_mkdir(transformations_dir);
        old_path = fullfile(source_patient_path,which_file);
        new_path = transformations_dir;
        coregDir = strrep(transformations_dir,'normalization/transformations','coregistration/anat');
        if exist(old_path,'file')
            %first move%
            copyfile(old_path,new_path);
            %then rename%
            if endsWith(which_file,'.h5') %already renames
                bids_name = [patient_name,'_',bids_name];
                if contains(which_file,'glanatComposite.h5')
                    reference = fullfile(ea_space,'t1.nii');
                elseif contains(which_file,'glanatInverseComposite.h5')
                    coregfiles = dir(fullfile(coregDir,'sub-*_ses-preop*'));
                    if isempty(coregfiles)
                       coregfiles = dir(fullfile(old_path,'*anat_*'));

                    end
                    reference = fullfile(coregfiles(1).folder,coregfiles(1).name);
                end
                try
                    ea_conv_antswarps(fullfile(new_path,which_file), reference);
                    outfile = strrep(fullfile(new_path,which_file),'.h5','.nii.gz');
                    movefile(outfile,fullfile(new_path,bids_name));
                catch
                    movefile(outfile,fullfile(new_path,bids_name));
                    warning('Transform files could not be converted to nii.gz format. Please check that the files are not corrupt manually.');
                end
            else
                disp(['Renaming file ' which_file ' to ' bids_name]);
                rename_path = fullfile(new_path,which_file);
                derivatives_cell{end+1,1} = fullfile(old_path);
                derivatives_cell{end,2} = fullfile(new_path,[patient_name,'_',bids_name]);
                movefile(rename_path,fullfile(new_path,[patient_name,'_',bids_name]));
            end
        end
    end


function move_raw2bids(source_patient_path,new_path,which_file,bids_name)
    tmp_path = fullfile(new_path,'tmp');
    ea_mkdir(tmp_path);
    if exist(fullfile(source_patient_path,which_file),'file')
        copyfile(fullfile(source_patient_path,which_file),tmp_path);
         if endsWith(which_file,'.nii')
            gzip(fullfile(tmp_path,which_file))
            ea_delete(fullfile(tmp_path,which_file))
            which_file = [which_file,'.gz'];
         end
        movefile(fullfile(tmp_path,which_file),fullfile(tmp_path,bids_name));
        copyfile(fullfile(tmp_path,bids_name),new_path);

    end
function move_mni2bids(mni_files,native_files,~,headmodel,which_pipeline,patient_name,new_path)
    if strcmp(which_pipeline,'headmodel')
        if ~isempty(mni_files)
            for headmodel_mni_file = 1:length(mni_files)
                if ismember(mni_files{headmodel_mni_file},headmodel{:,1})
                    indx = cellfun(@(x)strcmp(x,mni_files{headmodel_mni_file}),headmodel{:,1});
                    movefile(fullfile(new_path,which_pipeline,'MNI152NLin2009bAsym',mni_files{headmodel_mni_file}),fullfile(new_path,which_pipeline,'MNI152NLin2009bAsym',[patient_name,'_',headmodel{1,2}{indx}]));
                end
            end
        end
        if ~isempty(native_files)
            for headmodel_native_file = 1:length(native_files)
                if ismember(native_files{headmodel_native_file},headmodel{:,1})
                    indx = cellfun(@(x)strcmp(x,native_files{headmodel_native_file}),headmodel{:,1});
                    movefile(fullfile(new_path,which_pipeline,'native',native_files{headmodel_native_file}),fullfile(new_path,which_pipeline,'native',[patient_name,'_',headmodel{1,2}{indx}]));
                end
            end
        end

    end

function generate_rawImagejson(patient_name,dest)
    %prev inputs: doDcmConv,doOnlyRaw
    output_dir = fullfile(dest,'derivatives','leaddbs',patient_name,'prefs');
    ea_mkdir(output_dir);

    raw_data_path = fullfile(dest,'rawdata',patient_name);

    coreg_dir = fullfile(dest,'derivatives','leaddbs',patient_name,'coregistration','anat');
    raw_preop_dir = fullfile(raw_data_path,'ses-preop','anat');
    raw_postop_dir = fullfile(raw_data_path,'ses-postop','anat');
    preprocessing_dir = fullfile(dest,'derivatives','leaddbs',patient_name,'preprocessing','anat');
    ea_mkdir(coreg_dir);
    ea_mkdir(raw_preop_dir);
    ea_mkdir(raw_postop_dir);
    ea_mkdir(preprocessing_dir);

    opt.FileName = fullfile(dest,'derivatives','leaddbs',patient_name,'prefs',[patient_name,'_','desc-rawimages.json']);
    %special_case

    preop_files = dir_without_dots(fullfile(raw_preop_dir,'*.nii.gz'));
    preop_files = {preop_files.name};
    preop_modalities = {};
    for comparing_files = 1:length(preop_files)
        preop_file = regexprep(preop_files{comparing_files},'\.nii(\.gz)?$','');
        preop_mod = strsplit(preop_file,'_');
        preop_mod = preop_mod{end};
        preop_modalities{end+1} = preop_mod;
    end
    coreg_t1 = dir(fullfile(coreg_dir,'sub-*_ses-preop_space-anchorNative_*_acq-*_T1w.nii'));
    rawdata_t1 = dir(fullfile(raw_preop_dir,'sub-*_ses-preop_acq-*_T1w.nii'));
    if ~isempty(coreg_t1) && isempty(rawdata_t1)
        transfer_t1 = 1;
    else
        transfer_t1 = 0;
    end
    if isempty(preop_files) || transfer_t1
        %other preop files
        coreg_preop_files = dir(fullfile(coreg_dir,'sub-*_ses-preop_space-anchorNative_*_acq-*.nii'));
        coreg_preop_files = {coreg_preop_files.name};
        coreg_preop_files(endsWith(coreg_preop_files, '_mask.nii')) = [];
        for coreg_files = 1:length(coreg_preop_files)
            coreg_file = regexprep(coreg_preop_files{coreg_files}, '\.nii(\.gz)?$', '');
            coreg_mod = strsplit(coreg_file,'_');
            coreg_mod = coreg_mod{end};
            if ~ismember(coreg_mod,preop_modalities) || isempty(preop_modalities)
                coreg_preproc_name = [strrep(coreg_file,'space-anchorNative_',''),'.nii'];
                coreg_raw_preop = [strrep(coreg_file,'space-anchorNative_desc-preproc_',''),'.nii'];
                %change the filename
                %also copy this file to the preproc folder
                copyfile(fullfile(coreg_dir,[coreg_file '.nii']),fullfile(preprocessing_dir,[coreg_file '.nii']));
                movefile(fullfile(preprocessing_dir,[coreg_file '.nii']),fullfile(preprocessing_dir,coreg_preproc_name));
                copyfile(fullfile(coreg_dir,[coreg_file,'.nii']),fullfile(raw_preop_dir,[coreg_file,'.nii']));
                movefile(fullfile(raw_preop_dir,[coreg_file,'.nii']),fullfile(raw_preop_dir,coreg_raw_preop));
                gzip(fullfile(raw_preop_dir,coreg_raw_preop));
                ea_delete(fullfile(raw_preop_dir,coreg_raw_preop));
                preop_files{end+1} = coreg_raw_preop;
                %preop_files{end+1} = coreg_file;

            end

        end
    end

    if ~isempty(preop_files)
        for i=1:length(preop_files)
            json_val = regexprep(preop_files{i},'\.nii(\.gz)?$','');
            if contains(preop_files{i},'acq-')
                temp_tag = strsplit(json_val,'-');
                modality_str = temp_tag{end};
                rawdata_fieldname = modality_str;
                anat_files_selected.preop.anat.(rawdata_fieldname) = json_val;
            else
                temp_tag = strsplit(json_val,'preop_');
                rawdata_fieldname = temp_tag{end};
                anat_files_selected.preop.anat.(rawdata_fieldname) = json_val;
            end
        end
    end

    postop_files = dir_without_dots(fullfile(raw_postop_dir,'*.nii.gz'));
    postop_files = {postop_files.name};
    postop_modalities = {};
    for comparing_files = 1:length(postop_files)
        postop_file = regexprep(postop_files{comparing_files},'\.nii(\.gz)?$','');
        postop_mod = strsplit(postop_file,'-');
        postop_mod = postop_mod{end};
        postop_modalities{end+1} = postop_mod;
    end
    if isempty(postop_files)
        coreg_postop_files = dir(fullfile(coreg_dir,'sub-*_ses-postop_space-anchorNative_desc-preproc_*.nii'));
        coreg_postop_files = {coreg_postop_files.name};
        for coreg_files = 1:length(coreg_postop_files)
            coreg_file = regexprep(coreg_postop_files{coreg_files}, '\.nii(\.gz)?$', '');
            coreg_mod = strsplit(coreg_file,'-');
            coreg_mod = coreg_mod{end};
            if ~ismember(coreg_mod,postop_modalities) || isempty(postop_modalities)
                coreg_postop_name = [strrep(coreg_file,'space-anchorNative_',''),'.nii'];
                coreg_raw_postop = [strrep(coreg_file,'space-anchorNative_desc-preproc_',''),'.nii'];
                %change the filename
                %also copy this file to the preproc folder
                copyfile(fullfile(coreg_dir,[coreg_file '.nii']),fullfile(preprocessing_dir,[coreg_file '.nii']));
                movefile(fullfile(preprocessing_dir,[coreg_file '.nii']),fullfile(preprocessing_dir,coreg_postop_name));
                copyfile(fullfile(coreg_dir,[coreg_file,'.nii']),fullfile(raw_postop_dir,[coreg_file,'.nii']));
                movefile(fullfile(raw_postop_dir,[coreg_file,'.nii']),fullfile(raw_postop_dir,coreg_raw_postop));
                gzip(fullfile(raw_postop_dir,coreg_raw_postop));
                ea_delete(fullfile(raw_postop_dir,coreg_raw_postop));
                postop_files{end+1} = coreg_raw_postop;
                %postop_files{end+1} = coreg_file;
            end
        end

    end
    for i=1:length(postop_files)
        json_val = regexprep(postop_files{i},'\.nii(\.gz)?$','');
        if contains(postop_files{i},'ct','IgnoreCase',true)
            rawdata_fieldname = 'CT';
            anat_files_selected.postop.anat.(rawdata_fieldname) = json_val;
        else
            temp_tag = strsplit(json_val,'-');
            rawdata_fieldname = temp_tag{end};
            anat_files_selected.postop.anat.(rawdata_fieldname) = json_val;

        end
    end
    if exist('anat_files_selected','var')
        savejson('',anat_files_selected,opt);
    else
        ea_cprintf('CmdWinWarnings', 'No rawimages or coregistration files found, therefore no rawimages.json set!\n');
    end



function bids_name = add_tag(try_bids_name,mod_cell,tag_cell)
    if contains(try_bids_name, 'acqTag_mod-')
        bids_mod = regexp(try_bids_name, '(?<=_mod-)[a-zA-Z0-9]+', 'match', 'once');
    else
        bids_mod = regexp(try_bids_name, '(?<=acqTag_)[a-zA-Z0-9]+', 'match', 'once');
    end
    indx = cellfun(@(x)isequal(x,bids_mod),mod_cell);
    if isempty(find(indx,1))
        bids_name = strrep(try_bids_name,'acqTag_','');
    else
        tag = tag_cell{indx};
        bids_name = strrep(try_bids_name,'acqTag',['acq-' tag]);
    end


function bids_name = CheckifAlreadyExists(path,try_bids_name)
    if contains(try_bids_name,'acq-')
        split_bids_name = strsplit(try_bids_name,'acq-');
        tag_bids_name = split_bids_name{2};
        tag_bids_name = strsplit(tag_bids_name,'_');
        tag = tag_bids_name{1};
        bids_mod = tag_bids_name{2};
        if exist(fullfile(path,try_bids_name),'file')
            new_bids_name = [split_bids_name{1},'acq-',tag,num2str(1),'_',bids_mod];
            movefile(fullfile(path,try_bids_name),fullfile(path,new_bids_name));
            bids_name = [split_bids_name{1},'acq-',tag,num2str(2),'_',bids_mod];
        elseif exist(fullfile(path,[split_bids_name{1},'acq-',tag,num2str(1),'_',bids_mod]),'file')
            suffix = 2;
            filename_to_check = [split_bids_name{1},'acq-',tag,num2str(suffix),'_',bids_mod];
            while exist(fullfile(path,filename_to_check),'file')
                suffix = suffix + 1;
                filename_to_check = [split_bids_name{1},'acq-',tag,num2str(suffix),'_',bids_mod];
            end
            bids_name = filename_to_check;
        else
            bids_name = try_bids_name;
        end

    elseif endsWith(try_bids_name,'.mat')
        split_bids_name = strsplit(try_bids_name,'.mat');
        split_bids_name = split_bids_name{1};
        if exist(fullfile(path,try_bids_name),'file')
            new_bids_name = [split_bids_name,num2str(1),'.mat'];
            movefile(fullfile(path,try_bids_name),fullfile(path,new_bids_name));
            bids_name = [split_bids_name,num2str(2),'.mat'];
        elseif exist(fullfile(path,[split_bids_name,num2str(1),'.mat']),'file')
            suffix = 2;
            filename_to_check = [split_bids_name,num2str(suffix),'.mat'];
            while exist(fullfile(path,filename_to_check),'file')
                suffix = suffix + 1;
                filename_to_check = [split_bids_name,num2str(suffix),'.mat'];
            end
            bids_name = filename_to_check;
        else
            bids_name = try_bids_name;
        end
    else
        bids_name = try_bids_name;
    end


function [matching_files,noMatch] = match_exact(cell_files,expr)
    tmp_cell = cellfun(@(x) ~isempty(x) && x(1) == 1,regexpi(cell_files, expr));
    matching_files = cell_files(tmp_cell);
    noMatch = cell_files(~tmp_cell);


function bids_mod = add_mod(to_match,legacy_modalities,rawdata_containers)
    to_match = regexprep(to_match, '[\s]', '');
    for legacy_mod=1:length(legacy_modalities)
         modality_str = legacy_modalities{legacy_mod};
         if endsWith(to_match,'.nii') || endsWith(to_match,'.nii.gz')
             if contains(to_match,modality_str,'IgnoreCase',true)
                 bids_mod = rawdata_containers(modality_str);
                 break;
             elseif legacy_mod == length(legacy_modalities)
                 bids_mod = regexp(to_match, '(?<=^(raw_)?anat_).*(?=\.nii$)', 'match', 'once'); % file should have a [raw_]anat_ something catch
                 bids_mod = upper(bids_mod);
             end
         elseif endsWith(to_match,'.png')
             to_match_str = strsplit(to_match,'anat_t1');
             to_match_mod = to_match_str{1};
             if contains(to_match_mod,modality_str,'IgnoreCase',true)
                 bids_mod = rawdata_containers(modality_str);
                 break;
             elseif legacy_mod == length(legacy_modalities)
                 [~,to_match_mod,~] = fileparts(to_match_mod);
                 try
                     bids_mod = upper(to_match_mod(end-3:end));
                     bids_mod = regexprep(bids_mod, '[0-9_]', '');
                 catch
                     bids_mod = 'misc';
                 end
             end
         else
            bids_mod = '';
          end

    end


function model_name = add_model(stimFolder)
    stimParams = ea_regexpdir(stimFolder, 'stimparameters\.mat$', 0);
    if ~isempty(stimParams)
        load(stimParams{1},'S')
        model_name = ea_simModel2Label(S.model);
    else
        ea_cprintf('CmdWinWarnings', 'Missing stimparameters under %s\nSet to SimBio model by default, please check manually.\n', stimFolder);
        model_name = 'simbio';
    end


function generate_bidsConnectome_name(mni_folder,connectome_folder,lead_mapper,stimulations,tag_struct)
    %connectome_folder should be the full path of the connectome files,
    %mni_folder =
    %fullfile(new_path,pipeline,'MNI152NLinAsym','stimulations',stimulation
    %name)
    mapper_output_files = dir_without_dots(connectome_folder);
    mapper_output_files = {mapper_output_files(~[mapper_output_files.isdir]).name};
    ea_mkdir(mni_folder);

    for mapper_file = 1:length(mapper_output_files)
        matching_stim_file = regexprep(mapper_output_files{mapper_file},'(fl_|efield_|efield_gauss_)', '');
        matching_lm_file = ['_' matching_stim_file];
        if ismember(matching_lm_file,lead_mapper{:,1})
            if contains(mapper_output_files{mapper_file},'efield')
                tag_struct.simtag = 'efield';
            elseif contains(mapper_output_files{mapper_file},'efield_gauss')
                tag_struct.simtag = 'efieldgauss';
            else %doesn't contain either
                tag_struct.simtag = 'binary';
            end
            extract_old_hemisdesc_str = strsplit(mapper_output_files{mapper_file},'_');

            indx = cellfun(@(x)strcmp(x,matching_lm_file),lead_mapper{:,1});
            bids_name = lead_mapper{1,2}{indx};
            %replace hemidesc tag
            if strcmp(extract_old_hemisdesc_str{1},'fl')
                bids_name = strrep(bids_name,'flippedTag','flipped');
            else
                bids_name = strrep(bids_name,'_hemidesc-flippedTag','');
            end
            %replace sim tag
            bids_name = strrep(bids_name,'simTag',tag_struct.simtag);
            %replace model tag
            bids_name = strrep(bids_name,'modelTag',tag_struct.modeltag);
            %replace connectome tag
            bids_name = strrep(bids_name,'conTag',tag_struct.conntag);
            copyfile(fullfile(connectome_folder,mapper_output_files{mapper_file}),fullfile(mni_folder,mapper_output_files{mapper_file}));
            movefile(fullfile(mni_folder,mapper_output_files{mapper_file}),fullfile(mni_folder,[tag_struct.subjID,'_',bids_name]));
            fprintf('Renaming file %s as %s.\n',mapper_output_files{mapper_file},[tag_struct.subjID,'_',bids_name]);

        elseif ismember(matching_stim_file,stimulations{:,1})
             if contains(mapper_output_files{mapper_file},'efield')
                tag_struct.simtag = 'efield';
            elseif contains(mapper_output_files{mapper_file},'efield_gauss')
                tag_struct.simtag = 'efieldgauss';
            else %doesn't contain either
                tag_struct.simtag = 'binary';
            end
            extract_old_hemisdesc_str = strsplit(mapper_output_files{mapper_file},'_');

            indx = cellfun(@(x)strcmp(x,matching_stim_file),stimulations{:,1});
            bids_name = stimulations{1,2}{indx};
            %replace hemidesc tag
            if strcmp(extract_old_hemisdesc_str{1},'fl')
                bids_name = strrep(bids_name,'flippedTag','flipped');
            else
                bids_name = strrep(bids_name,'_hemidesc-flippedTag','');
            end
            %replace sim tag
            bids_name = strrep(bids_name,'simTag',tag_struct.simtag);
            %replace model tag
            bids_name = strrep(bids_name,'modelTag',tag_struct.modeltag);
            %replace connectome tag
            bids_name = strrep(bids_name,'conTag',tag_struct.conntag);
            copyfile(fullfile(connectome_folder,mapper_output_files{mapper_file}),fullfile(mni_folder,mapper_output_files{mapper_file}));
            movefile(fullfile(mni_folder,mapper_output_files{mapper_file}),fullfile(mni_folder,[tag_struct.subjID,'_',bids_name]));
            fprintf('Renaming file %s as %s.\n',mapper_output_files{mapper_file},[tag_struct.subjID,'_',bids_name]);
        else
            evalin('base','WARNINGSILENT=1;');
            ea_warning(sprintf('BIDS tag could not be assigned for %s. Please rename manually',mapper_output_files{mapper_file}));
        end

    end
    %evalin('base','WARNINGSILENT=1;');
    %ea_warning(sprintf('Deleting old copy of connectome folder %s. You can find it in the source patient folder if you need.',connectome_folder));
    %ea_delete(fullfile(connectome_folder));


function files_to_move = reorderfiles(files_to_move)
newfiles={};
k=1;
%send coregmethod and normmethod back
while k<= length(files_to_move)
    if endsWith(files_to_move{k},'.mat')
        newfiles{end+1} = files_to_move{k};
        files_to_move(k) = [];
    else
        k = k+1;
    end
end
files_to_move = [files_to_move; newfiles'];


function file2json(fname_in,fname_out,derivatives_cell)
legacy_modalities = {'t1','t2star','pd','ct','tra','cor','sag','fgatir','fa','dti','dti.bval','dti.bvec','t2','flair','inv','swi'};
bids_modalities = {'T1w','T2starw','PDw','CT','acq-ax_MRI','acq-cor_MRI','acq-sag_MRI','FGATIR','fa','dwi','dwi.bval','dwi.bvec','T2w','FLAIR','INV','SWI'};
rawdata_containers = containers.Map(legacy_modalities,bids_modalities);
opt.FileName = fname_out;
[filepath,filename,~] = fileparts(fname_in);
[op_dir,pt_name,~] = fileparts(fname_out);
pt_name = strsplit(pt_name,'_desc-');
pt_name = pt_name{1};
op_dir = strrep(op_dir,'log','anat');
json_mat = struct();
old = {derivatives_cell{:,1}};
new = {derivatives_cell{:,2}};
%function to convert mat files and text
if endsWith(fname_in,'.mat')
    %special case for legacy dataset
    %read the input mat
    %dealing with coregistration
    if strcmp(filename,'ea_coreg_approved')
        input_mat = load(fname_in);
        coreg_fieldnames = fieldnames(input_mat);
        %determine files inside the coreg files
        for i=1:length(coreg_fieldnames)
            if startsWith(coreg_fieldnames{i},'anat') || ~isempty(regexp(coreg_fieldnames{i}, '_tra|_sag|_cor', 'once'))
                if ismember(fullfile(filepath,[coreg_fieldnames{i},'.nii']),old)
                    indx = cellfun(@(x)strcmp(x,fullfile(filepath,[coreg_fieldnames{i},'.nii'])),old);
                    try_bids_name = new{indx};
                    [~,bids_name,~] = fileparts(try_bids_name);
                else
                    %try to convert to bids mod
                    bids_mod = add_mod([coreg_fieldnames{i},'.nii'],legacy_modalities,rawdata_containers);
                    if contains(bids_mod,'_MRI')
                        bids_name = [pt_name,'_ses-postop_space-anchorNative_desc-preproc_',bids_mod,'.nii'];
                    else
                        try
                            tag = ea_checkacq(fullfile(filepath,[coreg_fieldnames{i},'.nii']));
                        catch
                            warning("tag may not be set correctly, please recheck if you have issues with coregistration!");
                            tag = 'iso';
                        end
                        bids_name = [pt_name,'_ses-preop_space-anchorNative_desc-preproc_acq-',tag,'_',bids_mod,'.nii'];
                    end
                    bids_name = CheckifAlreadyExists(op_dir,bids_name);
                end
                field_name = strsplit(bids_name,'acq-');
                mod = regexprep(field_name{end},'\.nii$','');
                json_mat.approval.(mod) = input_mat.(coreg_fieldnames{i});
            elseif endsWith(coreg_fieldnames{i}, 'postop_ct')
                json_mat.approval.('CT') = input_mat.(coreg_fieldnames{i});
            end
        end
        if isfield(json_mat,'approval')
            savejson('',json_mat,'approval',json_mat.approval,opt);
        else
            warning("Coregistration files were not transformed. Please review your BIDSified folders with caution")
        end
    elseif strcmp(filename, 'ea_coregctmethod_applied')
        temp_mat = load(fullfile(filepath,'ea_coregctmethod_applied.mat'));
        method_used = generateMethod(temp_mat,'coregct_method_applied');
        if isempty(method_used)
            % Fallback to default method, should never reach here
            method_used = 'ANTs (Avants 2008)';
            ea_cprintf('CmdWinWarnings', 'CT coregistration method fallbacks to default: ''%s''.\n', method_used);
        end
        modality = 'CT';

        if isfile(fname_out)
            json_mat = loadjson(fname_out);
        end

        % Set approval status to 1 in case field is missing
        if ~isfield(json_mat, 'approval') || ~isfield(json_mat.approval, modality)
            json_mat.approval.(modality) = 1;
        end
        json_mat.method.(modality) = method_used;

        if isfield(json_mat,'method')
            savejson('',json_mat,'method',json_mat.method,opt);
        end
    elseif strcmp(filename,'ea_coregmrmethod_applied')
        temp_mat = load(fullfile(filepath,'ea_coregmrmethod_applied.mat'));
        method_used = generateMethod(temp_mat,'coregmr_method_applied');
        if isempty(method_used)
            % Fallback to default method, since coregmr method is not
            % properly stored in classic version of LeadDBS
            method_used = 'SPM (Friston 2007)';
            ea_cprintf('CmdWinWarnings', 'MR coregistration method fallbacks to default: ''%s''.\n', method_used);
        end

        if isfile(fname_out)
            json_mat = loadjson(fname_out);
        end

        if isfield(json_mat, 'approval')
            % If ea_coreg_approved.mat exists and was properly migrated,
            % use MR modalities from the *_desc-coregmethod.json first.
            modalities = setdiff(fieldnames(json_mat.approval), 'CT');
        end

        if isempty(modalities)
            % In case modalities not found above, check the migrated MRIs
            modalities = checkModalities(fullfile(fileparts(fileparts(fname_out)), 'anat'));
        end

        for m=1:length(modalities)
            % Set approval status to 1 in case field is missing
            if ~isfield(json_mat, 'approval') || ~isfield(json_mat.approval, modalities{m})
                json_mat.approval.(modalities{m}) = 1;
            end
            json_mat.method.(modalities{m}) = method_used;
        end

        if isfield(json_mat,'method')
            savejson('',json_mat,'method',json_mat.method,opt);
        end
    elseif strcmp(filename,'ea_normmethod_applied')
        input_mat = load(fname_in);
        normalize_fieldnames = fieldnames(input_mat);
        for i=1:length(normalize_fieldnames)
            method_used = generateMethod(input_mat,'norm_method_applied');
            json_mat.approval = 1;
            json_mat.method = method_used;
        end
        savejson('',json_mat,opt);
    end
end


function method_used = generateMethod(input_mat,modality_field)
if iscell(modality_field)
    modality_field = modality_field{1};
end

if contains(modality_field, 'coregct')
    regType = 'CT coregistration';
elseif contains(modality_field, 'coregmr')
    regType = 'MR coregistration';
elseif contains(modality_field, 'norm')
    regType = 'normalization';
end

if isfield(input_mat,modality_field)
    if ischar(input_mat.(modality_field))
        input_mat.(modality_field) = {input_mat.(modality_field)};
    end
    if strcmp(input_mat.(modality_field){end},'ea_normalize_apply_normalization')
        j = length(input_mat.(modality_field));
        for i=1:length(input_mat.(modality_field))
            if ~strcmp(input_mat.(modality_field){j},'ea_normalize_apply_normalization')
                modField = input_mat.(modality_field){j};
                %need to get the first entry thats not
                %ea_normalize_apply_normalization
                break
            end
            j = j-1;
        end
    else
        modField = input_mat.(modality_field){end};
    end
    if strcmp(modField,'ANTs') || contains(modField,'_ants')
        method_used = ea_normalize_ants('promt');
    elseif strcmp(modField,'BRAINSFit')
        method_used = 'BRAINSFit (Johnson 2007)';
    elseif strcmp(modField,'FLIRT')
        method_used = 'FLIRT (Jenkinson 2001 & 2002)';
    elseif strcmp(modField,'BBR')
        method_used = 'FLIRT BBR (Greve and Fischl 2009)';
    elseif strcmp(modField,'Hybrid SPM & ANTs')
        method_used = 'Hybrid SPM & ANTs';
    elseif strcmp(modField,'Hybrid SPM & BRAINSFIT')
        method_used = 'Hybrid SPM & BRAINSFIT';
    elseif strcmp(modField,'SPM') || contains(modField,'_spm')
        method_used = 'SPM (Friston 2007)';
    elseif strcmp(modField,'Hybrid SPM & FLIRT')
        method_used = 'Hybrid SPM & FLIRT';
    elseif strcmp(modField,'Schoenecker 2009')
        method_used = 'Three-step affine normalization (ANTs; Schonecker 2009)';
    else
        method_used = '';
        ea_cprintf('CmdWinWarnings', "Failed to identify the %s method used. Please take a closer look manually.\n", regType);
    end
else
    method_used = '';
    ea_cprintf('CmdWinWarnings', "Failed to identify the %s method used. Please take a closer look manually.\n", regType);
end


function fileList = filterLegacyFiles(fileList)
    startsWithPattern = {'lpostop', 'lanat', 'ea_pm', 'rc', 'u_rc', 'j_rc', ...
        'v_rc', 'y_anat_', 'iy_anat_', 'F.', 'DICOMDIR', 'Thumbs', 'Icon'};
    fileList(startsWith(fileList, startsWithPattern, 'IgnoreCase', true)) = [];

    endsWithPattern = {'.ps'};
    fileList(endsWith(fileList, endsWithPattern, 'IgnoreCase', true)) = [];

    containsPattern = {'tmp', 'temp', 'copy', 'grid', 'seg8', 'ignore', 'backup'};
    fileList(contains(fileList, containsPattern, 'IgnoreCase', true)) = [];

    fileList(endsWith(fileList, '.h5') & ~startsWith(fileList, 'glanat')) = [];


function modalities = checkModalities(coregAnatFolder)
    % In case modalities not found above, check the migrated MRIs
    coregAnat = ea_regexpdir(coregAnatFolder, 'acq-.*\.nii$', 0, 'f');
    modalities = unique(regexp(coregAnat, '(?<=_acq-).*(?=\.nii$)', 'match', 'once'));
    % Remove anchor modality
    prefs = ea_prefs;
    anchorModality = prefs.prenii_order{1};
    if sum(contains(modalities, anchorModality)) == 1
        % Only one anchor modality image exists
        modalities(contains(modalities, anchorModality)) = [];
    elseif sum(contains(modalities, anchorModality)) > 1
        % Multiple anchor modality images exist, remove the first
        % one according to the pre-defined order: iso, ax, cor, sag
        otherModalities = modalities(~contains(modalities, anchorModality));
        anchorModalities = setdiff(modalities, otherModalities);
        [~, ind] = ea_sortalike(lower(regexp(anchorModalities, '[^\W_]+(?=_[^\W_]+)', 'match', 'once')), {'iso', 'ax', 'cor', 'sag'});
        anchorModalities = anchorModalities(ind);
        modalities = [anchorModalities(2:end); otherModalities];
    end
