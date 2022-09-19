function ea_lead_import(source_folder,options,handles,dest_folder)


sub_counter = 1;
for pt = 1:length(source_folder)
    ea_busyaction('on',handles.leadfigure,'dbs');
    
    if exist('source_folder','var')
        count = checkValidity(source_folder{pt});
        if count >= 1
            if exist('handles','var')
                setappdata(handles.leadfigure, 'BIDSRoot', {});
                setappdata(handles.leadfigure, 'subjID', {});
            end
            if length(source_folder) == 1
                return
            else
                continue
            end
        end

        %detect dataset type
        [results,flag,dicom_conv,doMigrate,doOnlyRaw,dicom_source_folder] = enablerules(source_folder{pt},options);

        if strcmp(results,'fail')
            if length(source_folder) == 1
                if exist('handles','var')
                    setappdata(handles.leadfigure, 'BIDSRoot', {});
                    setappdata(handles.leadfigure, 'subjID', {});
                end
                return
            else
                continue
            end
        else
            [~,subj,~] = fileparts(source_folder{pt});
            subjID{sub_counter} =subj;
            sub_counter = sub_counter+1;
        end
    end
    if iscell(dest_folder)
        dest_folder = dest_folder{1};
    end
    if exist('handles','var')
        if ~isempty(getappdata(handles.leadfigure, 'BIDSRoot'))
            if contains(source_folder{1},'sourcedata')
                dest_folder = {getappdata(handles.leadfigure, 'BIDSRoot')};
            end
        end
    end
    run_prop.source = source_folder{pt};
    run_prop.dicom_source = dicom_source_folder;
    run_prop.dest = dest_folder;
    run_prop.convert_using = options.prefs.migrate.DicomConversionTool;
    run_prop.flag = flag;
    run_prop.dicom_bool = dicom_conv;
    run_prop.migrate_bool = doMigrate;
    run_prop.doOnlyRaw = doOnlyRaw;
    runImport(run_prop)
    
end
setappdata(handles.leadfigure, 'BIDSRoot', dest_folder);
setappdata(handles.leadfigure, 'subjID', subjID');
ea_busyaction('off',handles.leadfigure,'dbs');

function runImport(run_prop)
    source_dir = run_prop.source;
    dicom_source_folder = run_prop.dicom_source;
    dest_folder = run_prop.dest;
    convert_using = run_prop.convert_using;
    flag = run_prop.flag;
    dicom_conv = run_prop.dicom_bool;
    doMigrate = run_prop.migrate_bool;
    doOnlyRaw = run_prop.doOnlyRaw;
    if strcmp(flag,'onlyMigrate')
        disp('Migrating dataset...');
        if dicom_conv && ~isempty(dicom_source_folder) && ~doOnlyRaw
            ea_legacy2bids(source_dir,dest_folder,1,dicom_source_folder,1,0);
        elseif ~dicom_conv && doMigrate && doOnlyRaw
            ea_legacy2bids(source_dir,dest_folder,0,'',0,1);
        elseif doMigrate && ~dicom_conv && ~doOnlyRaw
            ea_legacy2bids(source_dir,dest_folder,0,'',0,0);
        end
    end
    if strcmp(flag,'MigrateDCMconv')
        disp('Migrating dataset...')
        %app.dicom_source: parent dicom,
        %implement checkpoint if folder is not empty
        if doOnlyRaw && dicom_conv
            ea_legacy2bids(source_dir,dest_folder,1,dicom_source_folder,1,1);
        elseif dicom_conv && ~doOnlyRaw
            ea_legacy2bids(source_dir,dest_folder,1,dicom_source_folder,1,0)
        end
        ea_dataset_import(source_dir,dest_folder,convert_using, 1);
    end
    % only DICOM conversion, no migration
    if strcmp(flag,'onlyDCMconv')
        ea_dataset_import(source_dir,dest_folder,convert_using, 1);
        
    end
    % if running automatically, close as soon as the migration
    % is over.

function count = checkValidity(folders)

    count=0;
   % for i=1:length(folders)
   path=folders;
   if isempty(path)
       return;
   end
   found = ea_checkSpecialChars(path);
   [~,dataset_folder,~] = fileparts(path);
   is_underscore = regexp(dataset_folder,'.*_.*','match');
   if ~isempty(is_underscore)
       ea_warndlg("It looks like you have underscores in your dataset folder name: %s. Please remove it and try again...",dataset_folder)
   end
   if found || ~isempty(is_underscore)
       count=count+1;
   end
   % end


function [results,flag,dicom_conv,doMigrate,doOnlyRaw,dicom_source_folder] = enablerules(selection,options)
    [filepath,subject_name,~] = fileparts(selection);
    
    subfolder = dir_without_dots(selection);
    subfolder = {subfolder.name};
    results = 'pass';
    dicom_source_folder = {};
    % detection function: if the patient has derivatives and
    % raw data, but not in BIDS format then it should be migrated.
    if endsWith(filepath, 'sourcedata')
        % Input is sourcedata folder under a bids dataset
        dicom_conv = 1;
        doMigrate = 0;
        doOnlyRaw = 0;
        dicom_source_folder{1} = selection;
        flag = 'onlyDCMconv';
    
        % BIDS Compliant raw data nifti files are available! in this case,
        % dcm -> bids conversion should require special handling.
        % Only the conversion from nifti -> lead bids will be performed
    elseif checkIfOneExist(subfolder,'^glanat.*.nii$') && checkIfOneExist(subfolder,'anat.*.nii') && checkIfNoneExist(subfolder,'DICOM') && checkIfNoneExist(subfolder,'.*.dcm')
        dicom_conv = 0;
        doMigrate = 1;
        doOnlyRaw = 0;
        flag = 'onlyMigrate';
        % only raw nifti available, no dicom is available. Migrate
        % will be done, but only to move & rename the raw files
    elseif checkIfOnlyExist(subfolder,'.*.nii||.*.nii.gz') && ~any(ismember(subfolder,'ea_ui.mat'))
        dicom_conv = 0;
        doMigrate = 1;
        if any(ismember(subfolder,'rpostop_ct.nii')) || any(ismember(subfolder,'tp_rpostop_ct.nii')) || any(ismember(subfolder,'tp_glpostop_ct'))
            doOnlyRaw = 0;
        else
            doOnlyRaw = 1;
        end
        flag = 'onlyMigrate';
        % both dicom & derivatives available. Requires special
        % handling in the migrate code.
    elseif any(ismember(subfolder,'ea_ui.mat')) && checkIfOneExist(subfolder,'^glanat*.nii$') && checkIfOneExist(subfolder,'dicom||.*.dcm')
        doMigrate = 1;
        dicom_conv = 1;
        doOnlyRaw = 0;
        flag = 'MigrateDCMconv';
        if isfolder(fullfile(selection,'DICOM'))
            dicom_source_folder{1} = fullfile(selection,'DICOM');
        else
            if ~exist(fullfile(selection,'DICOM'),'dir')
                mkdir(fullfile(selection,'DICOM'))
            end
            for i=1:length(subfolder)
                if regexp(subfolder{i},'.*.dcm')
                    movefile(fullfile(selection,subfolder{i}),fullfile(selection,'DICOM'))
                end
            end
        end
        dicom_source_folder{1} = fullfile(selection,'DICOM');
    elseif any(ismember(subfolder,'ea_ui.mat')) && checkIfOneExist(subfolder,'[^anat]*.nii$')
        doMigrate = 1;
        dicom_conv = 0;
        doOnlyRaw = 1;
        flag = 'onlyMigrate';
    
    elseif checkIfOneExist(subfolder,'.*.nii||.*.nii.gz') && checkIfOneExist(subfolder,'dicom||.*.dcm') % raw nifti + dicom
        doMigrate = 1;
        dicom_conv = 1;
        doOnlyRaw = 1;
        flag = 'MigrateDCMconv';
        if isfolder(fullfile(selection,'DICOM'))
            dicom_source_folder{1} = fullfile(selection,'DICOM');
        else
            if ~exist(fullfile(selection,'DICOM'),'dir')
                mkdir(fullfile(selection,'DICOM'))
            end
            for i=1:length(subfolder)
                if regexp(subfolder{i},'.*.dcm')
                    movefile(fullfile(selection,subfolder{i}),fullfile(selection,'DICOM'))
                end
            end
            dicom_source_folder{1} = fullfile(selection,'DICOM');
        end
    
        % only DICOM files available
    elseif checkIfOneExist(subfolder,'dicom*||.*.dcm')
        %todo: remove support for .dcm files inside
        flag = 'onlyDCMconv';
        doMigrate = 0;
        dicom_conv = 1;
        doOnlyRaw = 0;
        if ~isempty(ea_regexpdir(selection, 'DICOM', 0, 'd'))
            dcm_folder = ea_regexpdir(selection, 'DICOM', 0, 'd');
            [~,dicom_folder,~] = fileparts(dcm_folder);
            if ~strcmp(dicom_folder,'DICOM')
    
                dicom_conv = 1;
                doOnlyRaw = 0;
                folders_inside = dir_without_dots(dcm_folder);
                subfolders_inside = folders_inside([folders_inside.isdir]); %only get the folders, since they have to be moved to DICOM]
    
                if ~exist(fullfile(selection,'DICOM'),'dir')
                    mkdir(fullfile(selection,'DICOM'))
                end
                subfolder_names = {subfolders_inside.name};
                for i=1:length(subfolder_names)
                    movefile(fullfile(selection,dicom_folder,subfolder_names{i}), fullfile(selection,'DICOM',subfolder_names{i}))
                end
                dicom_source_folder{1} = fullfile(selection,'DICOM');
    
            else
                if ~exist(fullfile(selection,'DICOM'),'dir')
                    mkdir(fullfile(selection,'DICOM'))
                end
                for i=1:length(subfolder)
                    if regexp(subfolder{i},'.*.dcm')
                        movefile(fullfile(selection,subfolder{i}),fullfile(selection,'DICOM'))
                    end
                end
                dicom_source_folder{1} = fullfile(selection,'DICOM');
            end
        end
    
    elseif checkIfOneExist(subfolder,'^ea_.*.mat$') || ea_ismember(subfolder,"current_headmodel") || ea_ismember(subfolder,"stimulation") % for detached files only, experimental
        dicom_conv = 0;
        doMigrate = 1;
        doOnlyRaw = 0;
        flag = 'onlyMigrate';
        % additional check for subfolders of dicom. Not just
        % folders inside which there might be dicom files.
    elseif ~isfolder(fullfile(selection,'DICOM')) && ~isfolder(fullfile(selection,'dicom'))
        is_dicom_inside = dir(fullfile(selection,'**/*.dcm'));
        if ~isempty(is_dicom_inside)
            dicom_conv = 1;
            doOnlyRaw = 0;
            folders_inside = dir_without_dots(selection);
            subfolders_inside = folders_inside([folders_inside.isdir]); %only get the folders, since they have to be moved to DICOM]
            if ~exist(fullfile(selection,'DICOM'),'dir')
                mkdir(fullfile(selection,'DICOM'))
            end
            subfolder_names = {subfolders_inside.name};
            for i=1:length(subfolder_names)
                movefile(fullfile(selection,subfolder_names{i}), fullfile(selection,'DICOM',subfolder_names{i}))
            end
            dicom_source_folder{1} = fullfile(selection,'DICOM');
            flag = 'onlyDCMconv';
        else
            ea_warndlg('%s: unable to detect the type of processing required, skipping for now', subject_name{1});
            results = 'skip';
        end
    
    else
        results = 'fail';%case not handled
    end
    %compare the results with the user defined prefs:
    if  ~options.prefs.migrate.doDicomConversion
        if dicom_conv
            if ~strcmp(flag,'onlyDCMconv')
                dicom_conv = 0;
                if strcmp(flag,'MigrateDCMconv')
                    flag = 'onlyMigrate';
                end
            else
                results = 'fail';
                warning("Only dicom files found for this patient, but your preference suggests to ignore dicom conversion. Skipping for now...reset prefs to include this patient");
            end
        end
    end



function result = checkIfOneExist(cell_of_filenames,str_to_match)
    result = 0;
    %regexp returns an array of cells ({0} = does not match str,
    %{[1]} = matches the str).
    %isempty returns an array of 1 when the regexp is not
    %satisfied. Inverse of the cell fun will return 1 if the regexp
    %matches the str, and 0 if regexp does not match the str
    %any will check if there is atleast one pattern has matched
    if any(~cellfun('isempty',regexpi(cell_of_filenames,str_to_match)))
        result = 1;
    end

function result = checkIfNoneExist(cell_of_filenames,str_to_match)
    output = checkIfOneExist(cell_of_filenames,str_to_match);
    result = ~output;
        
function result = checkIfOnlyExist(cell_of_filenames,str_to_match)
    result = 0;
    if all(~cellfun('isempty',regexpi(cell_of_filenames,str_to_match)))
        result = 1;
    end
function is_member_flag = ea_ismember(cellA,stringB)
    %%this function returns a boolean value instead of a cell
    is_member_flag = 0;
    if ~iscell(cellA)
        ea_error("You have to supply a cell as the first argument");
    end
    if ~isstring(stringB)
        ea_error("You have to supply a string as the second argument");
    end
    cell_output = ismember(cellA,stringB);
    boolen_value = find(cell_output);
    if length(boolen_value) >= 1
        is_member_flag = 1;
    end

return

