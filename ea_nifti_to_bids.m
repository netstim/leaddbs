function [sortedFiles, returnCode] = ea_nifti_to_bids(niiFiles, dataset_folder, subjID, preset)

% Function creates a GUI in order to select which files should be used to create a BIDS compliant dataset and
% which files should be used by lead-dbs. After selection, a rawdata folder will be created and selected files will be
% copy/pasted to the appropriate locations
%
% input
%   niiFiles (string/cell): folder/files to be loaded (files should be in the same directory)
%   dataset_folder (string): BIDS dataset folder
%   subjID (string): subjID, including the 'sub' tag
%
% output
%   sortedFiles (struct): struct with preop and postop session fields, inside each session field there are fields for every modality
%                                 those fields have strings with the filenames that are selected and to be used by lead-dbs
%   returnCode: 'okay', 'cancel' or 'close'
% __________________________________________________________________________________
% Copyright (C) 2021 Charite University Medicine Berlin, Movement Disorders Unit
% Johannes Achtzehn
% Part of this code is inspired by the dicm2nii tool of xiangrui.li@gmail.com

%% preparations
% lookup table to pre-allocate options for the table.% check whether pref exists, if not create it
if ~ispref('dcm2bids', 'lookuptable')
    lookup_table = loadjson(fullfile(ea_getearoot(), 'helpers', 'dicom_bids_lookuptable.json'));
    addpref('dcm2bids', 'lookuptable', lookup_table)
else
    lookup_table = getpref('dcm2bids', 'lookuptable');
end

% nii folder instead of nii files provided
if isfolder(niiFiles)
    niiFiles = ea_regexpdir(niiFiles, '\.nii(\.gz)$', 0);
    if isempty(niiFiles)
        ea_cprintf('CmdWinWarnings', 'NIfTI not found in the specified folder!\n')
        return;
    end
end

% gzip *.nii
for i=1:length(niiFiles)
    if endsWith(niiFiles{i}, '.nii')
        gzip(niiFiles{i});
        delete(niiFiles{i});
        niiFiles{i} = [niiFiles{i}, '.gz'];
    end
end

jsonFiles = regexprep(niiFiles, '\.nii(\.gz)?$', '.json');

h_wait = waitbar(0, 'Please wait while NIfTI images are being loaded');
read_nii = nii_viewer('func_handle', 'read_nii');
% cant initialize imgs and img_resolutions since some files might be corrupt /
% not readable so we dont know the definitive length. This is fine and does
% not lead to speed issues. / AH 2023
cnt=1;
for image_idx = 1:length(niiFiles)
    imgs{cnt,1} = struct();
    try
        [imgs{cnt,1}.p, imgs{cnt,1}.frm, imgs{cnt,1}.rg, imgs{cnt,1}.dim] = read_nii(niiFiles{cnt}, [], 0);

        % get .json and read it if possible
        if isfile(jsonFiles{cnt})
            try
                imgs{cnt,1}.json_sidecar = loadjson(jsonFiles{cnt});
                imgs{cnt,1}.json_found = 1;
            catch
                ea_cprintf('CmdWinWarnings', 'There was a problem while loading the .json file , please ensure correct .json format:\n%s\n', jsonFiles{cnt})
                imgs{cnt,1}.json_found = 0;
            end
        else
            imgs{cnt,1}.json_found = 0;
        end
    catch
        ea_cprintf('CmdWinWarnings', 'There was a problem while loading the .nii file, please ensure this is a correct .nii image:\n%s\n', niiFiles{cnt})

        % delete files that have failed to load from the list
        niiFiles(cnt) = [];       
        jsonFiles(cnt) = [];

        % also delete them from the preset
        if exist('preset', 'var')
            preset(cnt) = [];
        end

        continue; % some images (e.g. 1-dimensional) will not be readable, skip them.
    end
    imgs{cnt,1}.percentile = prctile(imgs{cnt,1}.p.nii.img(:), 95, 'all');
    imgs{cnt,1}.img_thresholded = imgs{cnt,1}.p.nii.img(:);
    imgs{cnt,1}.img_thresholded(imgs{cnt,1}.img_thresholded > imgs{cnt,1}.percentile) = nan;
    imgs{cnt,1}.img_thresholded(imgs{cnt,1}.img_thresholded < median(imgs{cnt,1}.p.nii.img(:))) = nan;
    imgs{cnt,1}.img_thresholded = imgs{cnt,1}.img_thresholded(~isnan(imgs{cnt,1}.img_thresholded));
    imgs_resolution{cnt,1} = [imgs{cnt,1}.p.pixdim(1), imgs{cnt,1}.p.pixdim(2), imgs{cnt,1}.p.pixdim(3)];

    cnt=cnt+1;
    waitbar(image_idx / length(niiFiles), h_wait, sprintf('Please wait while Niftii images are being loaded (%i/%i)', image_idx, length(niiFiles)));
end
close(h_wait);
% refresh json files only from the readable niftis

N_fnames = length(niiFiles);

if N_fnames == 0
    ea_cprintf('CmdWinWarnings', 'No compatible nifti images could be found. Please check if the images are 3-D image files.\n');
    sortedFiles = [];
    returnCode = 'cancel';
    return;
end

anat_modalities = {'T1w', 'T2w', 'FGATIR', 'FLAIR', 'T2starw', 'PDw'};  % a list of all supported modalities
func_dwi_modalities = {'bold', 'sbref', 'dwi'};
postop_modalities = {'CT', 'MRI'};          % specifically a list of modalities required for postoperative sessions, will be used to check if postop modalities have been found
postop_acq_tags = {'ax', 'cor', 'sag'};     % a list of required acq-tags for the postop MRI images

% options that should appear in the table
table_options = struct;
table_options.Session = {'preop', 'postop'};
table_options.Type = {'anat', 'func', 'dwi'};
table_options.Modality = [anat_modalities, func_dwi_modalities, postop_modalities];
table_options.Acquisition = [postop_acq_tags, 'iso'];

nifti_table = structfun(@(x) categorical(repmat({'-'}, [N_fnames,1]), ['-' x]), table_options, 'uni', 0);
if ~exist('preset', 'var')
    nifti_table.Include = false(N_fnames,1);
else
    nifti_table.Include = logical(preset);
end
[~, nifti_table.Filename] = cellfun(@ea_niifileparts, niiFiles, 'Uni', 0);
nifti_table.Run = repmat("-", [N_fnames,1]);
nifti_table.Task = repmat("-", [N_fnames,1]);

nifti_table = orderfields(nifti_table, [5 6 1 7 8 2 3 4]);
nifti_table = struct2table(nifti_table);

try
    nifti_table_preallocated = preallocate_table(nifti_table, lookup_table, imgs_resolution);
catch
    disp('Preallocation of table failed, defaulting to skip!');
    nifti_table_preallocated = nifti_table;
end

% create GUI
uiapp = nifti_to_bids;
setappdata(uiapp.UIFigure, 'niiFolder', fileparts(niiFiles{1}));

% populate table
uiapp.niiFileTable.Data = nifti_table_preallocated;
uiapp.niiFileTable.ColumnEditable = [true false true true true true true true];
uiapp.niiFileTable.UserData.task_set = zeros(N_fnames, 1);   % this will be used to only set task to 'rest' automatically one time (the first time), potentially other fields in the future
uiapp.niiFileTable.UserData.fnames = repmat("", [N_fnames,1]);   % this will be used to keep track of the filenames

% set subject ID and file path
uiapp.SubjectIDLabel.Text = sprintf('Conversion for subject: %s',subjID);
uiapp.FilepathLabel.Text = sprintf('BIDS root: %s',dataset_folder);

% update preview tree and expand it
uiapp.previewtree_subj.Text = subjID;
expand(uiapp.Tree, 'all');

preview_nii(uiapp, imgs, []); % set initial image to the first one
uiapp.niiFileTable.SelectionType = 'row';
uiapp.niiFileTable.Selection = 1;

%% set callbacks of main GUI
cell_change_callback(uiapp, subjID, imgs, anat_modalities, postop_modalities, postop_acq_tags, []) % call preview tree updater to get preallocated changes

uiapp.niiFileTable.CellSelectionCallback = @(src,event) preview_nii(uiapp, imgs, event); % callback for table selection -> display current selected image
uiapp.niiFileTable.CellEditCallback = @(src,event) cell_change_callback(uiapp, subjID, imgs, anat_modalities, postop_modalities, postop_acq_tags, event); % callback for cell change -> update uiapp tree and adjacent cells

uiapp.UIFigure.WindowScrollWheelFcn = @(src, event) scroll_nii(uiapp, event);     % callback for scrolling images

% OK button behaviour
uiapp.OKButton.ButtonPushedFcn = @(btn,event) ok_button_function(uiapp, table_options, dataset_folder, subjID, postop_modalities, postop_acq_tags);

% cancel button behaviour
uiapp.CancelButton.ButtonPushedFcn =  @(btn,event) cancel_button_function(uiapp);

% looup table behaviour
uiapp.LookupButton.ButtonPushedFcn = @(btn,event) lookup_button_function(uiapp, imgs, imgs_resolution, table_options, subjID, anat_modalities, postop_modalities);

setappdata(groot, 'sortedFiles', []);
setappdata(groot, 'returnCode', 'close');

waitfor(uiapp.UIFigure);

sortedFiles = getappdata(groot, 'sortedFiles');
returnCode = getappdata(groot, 'returnCode');

if ~isempty(sortedFiles)
    field_names = fieldnames(sortedFiles);
    empty_fields = cellfun(@(x) isempty(sortedFiles.(x)), field_names);
    remove_fields = field_names(empty_fields);
    sortedFiles = rmfield(sortedFiles, remove_fields);
end

end

%% define callback of man GUI
%% lookup button function
function lookup_button_function(uiapp, imgs, imgs_resolution, table_options, subjID, anat_modalities, postop_modalities)

lookup_table_gui = ea_default_lookup;
lookup_table = getpref('dcm2bids', 'lookuptable');

% populate table with current preferences
lookup_table_data = convert_lookup_struct_to_table(lookup_table);
lookup_table_gui.UITable.Data = lookup_table_data;

% load defaults button
lookup_table_gui.LoaddefaultsButton.ButtonPushedFcn = @(btn,event) load_default_lookups(lookup_table_gui);

% load .json button
lookup_table_gui.LoadjsonButton.ButtonPushedFcn = @(btn,event) load_json_file(lookup_table_gui);

% cancel button
lookup_table_gui.CancelButton.ButtonPushedFcn = @(btn,event) cancel_lookup_function(lookup_table_gui);

% save button
lookup_table_gui.SaveButton.ButtonPushedFcn = @(btn,event) save_lookup_function(uiapp, lookup_table_gui, imgs, imgs_resolution,  table_options, subjID, anat_modalities, postop_modalities);

waitfor(lookup_table_gui.UIFigure);
end

function save_lookup_function(main_gui, lookup_table_gui, imgs, imgs_resolution, table_options, subjID, anat_modalities, postop_modalities)

lookup_table = convert_table_to_lookup_struct(lookup_table_gui.UITable.Data);

setpref('dcm2bids', 'lookuptable', lookup_table);

T_preallocated = preallocate_table(main_gui.niiFileTable.Data, lookup_table, imgs_resolution);

main_gui.niiFileTable.Data = T_preallocated;

cell_change_callback(main_gui, subjID, imgs, anat_modalities, postop_modalities, postop_acq_tags, [])

delete(lookup_table_gui);

end

function table_struct = convert_table_to_lookup_struct(table_strings)

table_struct = struct();

for row_idx = 1:height(table_strings)

    keywords_string = strsplit(table_strings.Keywords(row_idx), ',');

    keywords_cell = cell(1, length(keywords_string));
    for word = 1:length(keywords_string)
        keywords_cell{1, word} = keywords_string(word);
    end

    table_struct.(table_strings.Type(row_idx)).(table_strings.Session(row_idx)).(table_strings.Modality(row_idx)) = keywords_cell;
end

end

function cancel_lookup_function(lookup_table_gui)
s = uiconfirm(lookup_table_gui.UIFigure, 'Do you really want to cancel loopup modification?', 'Confirm close', ...
    'Options', {'Yes', 'No'}, 'Icon', 'question');
switch s
    case 'Yes'
        delete(lookup_table_gui);
end
end

function load_json_file(lookup_table_gui)
[FileName,PathName,~] = uigetfile('.json');
fprintf('Loading lookup_table from .json file at %s\n', fullfile(PathName, FileName));
lookup_table = loadjson(fullfile(PathName, FileName));
lookup_table_data = convert_lookup_struct_to_table(lookup_table);
lookup_table_gui.UITable.Data = lookup_table_data;
end

function load_default_lookups(lookup_table_gui)
fprintf('Loading defaults from .json file at %s\n', fullfile(ea_getearoot(), 'helpers', 'dicom_bids_lookuptable.json'));
lookup_table = loadjson(fullfile(ea_getearoot(), 'helpers', 'dicom_bids_lookuptable.json'));
setpref('dcm2bids', 'lookuptable', lookup_table)
lookup_table_data = convert_lookup_struct_to_table(lookup_table);
lookup_table_gui.UITable.Data = lookup_table_data;
end

%% cell change callback
function cell_change_callback(uiapp, subjID, imgs, anat_modalities, postop_modalities, postop_acq_tags, event)

uiapp.previewtree_preop_anat.Children.delete;     % delete children
uiapp.previewtree_postop_anat.Children.delete;    % delete children

uiapp.previewtree_preop_dwi.Children.delete;      % delete children
uiapp.previewtree_postop_dwi.Children.delete;     % delete children

uiapp.previewtree_preop_func.Children.delete;     % delete children
uiapp.previewtree_postop_func.Children.delete;    % delete children

% go through all of the items and automatically set columns based on others
for i = 1:height(uiapp.niiFileTable.Data)
    modality = char(uiapp.niiFileTable.Data.Modality(i));

    if ~isempty(event)
        % check only for the current selected row and exclude "Filename" and "Include" column
        if event.Indices(1) ~= i || event.Indices(2) < 3
            continue;
        end

        % set modality for anat image automatically when switching session to "preop"
        if event.Indices(2) == 3 && event.NewData == "preop" && uiapp.niiFileTable.Data.Type(i) == "anat" && strcmp(modality, 'MRI')
            uiapp.niiFileTable.Data.Modality(i) = 'T1w';
            uialert(uiapp.UIFigure, 'Modality falls back to ''T1w'' by default when session changed to ''preop''. Please update the modality accordingly!', '');

        % set modality for anat image automatically when switching session to "postop"
        elseif event.Indices(2) == 3 && event.NewData == "postop" && uiapp.niiFileTable.Data.Type(i) == "anat" && ~strcmp(modality, 'CT')
            uiapp.niiFileTable.Data.Modality(i) = 'MRI';

        % set task and modality for anat image automatically when switching type to "anat"
        elseif event.Indices(2) == 6 && event.NewData == "anat"
            uiapp.niiFileTable.Data.Task(i) = '-';
            uiapp.niiFileTable.UserData.task_set(i) = 0;
            uialert(uiapp.UIFigure, 'Modality falls back to ''T1w'' by default when type changed to ''anat''. Please update the modality accordingly!', '');
            uiapp.niiFileTable.Data.Modality(i) = 'T1w';

        % set task and modality for func image automatically when switching type to "func"
        elseif event.Indices(2) == 6 && event.NewData == "func"
            uiapp.niiFileTable.Data.Task(i) = 'rest';
            uiapp.niiFileTable.UserData.task_set(i) = 1;
            uiapp.niiFileTable.Data.Modality(i) = 'bold';

        % set task and modality for dwi image automatically when switching type to "dwi"
        elseif event.Indices(2) == 6 && event.NewData == "dwi"
            uiapp.niiFileTable.Data.Task(i) = '-';
            uiapp.niiFileTable.UserData.task_set(i) = 0;
            uiapp.niiFileTable.Data.Modality(i) = 'dwi';

        % set type, session, run and task for preop anat image automatically
        elseif ismember(modality, anat_modalities)
            uiapp.niiFileTable.Data.Type(i) = 'anat';
            uiapp.niiFileTable.Data.Session(i) = 'preop';
            uiapp.niiFileTable.Data.Task(i) = '-';
            uiapp.niiFileTable.UserData.task_set(i) = 0;

        % set type, session, run and task for postop anat image automatically
        elseif ismember(modality, postop_modalities)
            uiapp.niiFileTable.Data.Type(i) = 'anat';
            uiapp.niiFileTable.Data.Session(i) = 'postop';
            uiapp.niiFileTable.Data.Task(i) = '-';
            uiapp.niiFileTable.UserData.task_set(i) = 0;

            % if current acquistion tag is not ax, cor or sag, reset it
            if ~ismember(uiapp.niiFileTable.Data.Acquisition(i), postop_acq_tags) && strcmp(modality, 'MRI')
                if uiapp.niiFileTable.Data.Acquisition(i) ~= "-"
                    uialert(uiapp.UIFigure, 'For postop MRIs, the acquisition tag may only be set to ''ax'' (must have) , ''sag'' or ''cor''.', 'Invalid acquisition tag for postop MRI');
                end
                uiapp.niiFileTable.Data.Acquisition(i) = 'ax';
            end

        % set type to func for bold modality
        elseif strcmp(modality, 'bold') && uiapp.niiFileTable.UserData.task_set(i) == 1
            uiapp.niiFileTable.Data.Type(i) = 'func';

        % set task to rest if bold is selected and task is not set and has not been automatically set before (do this only once)
        elseif strcmp(modality, 'bold') && uiapp.niiFileTable.UserData.task_set(i) == 0
            uiapp.niiFileTable.Data.Type(i) = 'func';
            uiapp.niiFileTable.Data.Task(i) = 'rest';
            uiapp.niiFileTable.UserData.task_set(i) = 1;

        % set type to dwi for dwi modality
        elseif strcmp(modality, 'dwi')
            uiapp.niiFileTable.Data.Type(i) = 'dwi';
            uiapp.niiFileTable.Data.Task(i) = '-';
        end
    else    % automatically set rest as task if func is detected (this is not covered by the lookupable as of yet)
        if strcmp(modality, 'bold') && uiapp.niiFileTable.UserData.task_set(i) == 0
            uiapp.niiFileTable.Data.Type(i) = 'func';
            uiapp.niiFileTable.Data.Task(i) = 'rest';
            uiapp.niiFileTable.UserData.task_set(i) = 1;
        end
    end
end

% go through all the ones that are not included but have session, modality and type set and enable them
for i = find(~uiapp.niiFileTable.Data.Include)'
    session = char(uiapp.niiFileTable.Data.Session(i));
    type = char(uiapp.niiFileTable.Data.Type(i));
    modality = char(uiapp.niiFileTable.Data.Modality(i));

    if ~isempty(event) % check this only for the current selected one
        if ~ismember('-', {session, type, modality}) && event.Indices(2) > 2 && event.Indices(1) == i
            uiapp.niiFileTable.Data.Include(i) = true;
        end
    end
end

% finally, go through all the files that have been selected to include and update them in the uitree
for i = find(uiapp.niiFileTable.Data.Include)'

    session = char(uiapp.niiFileTable.Data.Session(i));
    type = char(uiapp.niiFileTable.Data.Type(i));
    run = char(uiapp.niiFileTable.Data.Run(i));
    task = char(uiapp.niiFileTable.Data.Task(i));
    modality = char(uiapp.niiFileTable.Data.Modality(i));

    if strcmp(modality, "CT")
        uiapp.niiFileTable.Data.Acquisition(i) = "-";
    else
        % Ensure MRI has a proper acq tag set
        if uiapp.niiFileTable.Data.Acquisition(i) == "-"
            acq = ea_checkacq(imgs{i}.p.pixdim);
            if ~ismember(acq, postop_acq_tags) && strcmp(modality, 'MRI')
                ea_cprintf('CmdWinWarnings', '''%s'' acquisition tag detected.\nHowever, it may only be set to ''ax'' (must have) , ''sag'' or ''cor'' for postop MRIs.\nReset to ''ax'' now.\n', acq);
                acq = "ax";
            end
            uiapp.niiFileTable.Data.Acquisition(i) = acq;
        end
    end
    acq = char(uiapp.niiFileTable.Data.Acquisition(i));

    if ~ismember('-', {session, type, modality})  % check whether everything has been properly defined befor updating uitree

        fname = generate_bids_filename(subjID, session, run, task, acq, modality, []);
        ui_field = ['previewtree_' session '_' type];

        % now check if there are duplicate filenames in the tree
        if ~isempty(uiapp.(ui_field).Children) && any(ismember(fname, {uiapp.(ui_field).Children.Text}))

            % if duplicate file names have been found, also update the already present filename
            row_idx_duplicate_previewtree = find(cellfun(@(c) ischar(c) && strcmp(c, fname), {uiapp.(ui_field).Children.Text}));    % find the row of the other duplicate in the previewtree
            row_idx_duplicate_filetable = find_dupl_file(uiapp.niiFileTable, i, session, type, run, task, modality);

            % for bold and dwi data, try to sort out duplicates automatically with direction tag
            if strcmp(modality, 'bold') || strcmp(modality, 'dwi') && imgs{i, 1}.json_found == 1 && ~isempty(row_idx_duplicate_filetable)

                % get phase encoding for the current row
                ped_current_row = get_phase_encoding_direction(imgs{i, 1}.json_sidecar);
                ped_dupl_row = get_phase_encoding_direction(imgs{row_idx_duplicate_filetable, 1}.json_sidecar);

                if ~strcmp(ped_current_row, '') && ~strcmp(ped_dupl_row, '') && ~strcmp(ped_current_row, ped_dupl_row)   % if encoding directions have been found
                    fname = generate_bids_filename(subjID, session, run, task, acq, modality, ped_current_row);
                    uiapp.(ui_field).Children(row_idx_duplicate_previewtree).Text = generate_bids_filename(subjID, session, run, task, acq, modality, ped_dupl_row);
                    uiapp.niiFileTable.UserData.fnames(row_idx_duplicate_filetable) = generate_bids_filename(subjID, session, run, task, acq, modality, ped_dupl_row);
                else
                    fname = ['>> ', fname, ' <<'];  % change filename to indicate duplicate
                    uiapp.(ui_field).Children(row_idx_duplicate_previewtree).Text = fname;  % set the other duplicate to this filename as well
                    uiapp.niiFileTable.UserData.fnames(row_idx_duplicate_filetable) = fname;
                end
                % for others, just set the filename with >><< to warn the user
            else
                fname = ['>> ', fname, ' <<'];  % change filename to indicate duplicate
                uiapp.(ui_field).Children(row_idx_duplicate_previewtree).Text = fname;  % set the other duplicate to this filename as well
            end
        end

        uiapp.niiFileTable.UserData.fnames(i) = fname;
        uitreenode(uiapp.(ui_field), 'Text', fname);    % set the filename of the current row
    end

end

end


%% preallocate table on the left
function table_preallocated = preallocate_table(table, lookup_table, imgs_resolution)

table_preallocated = table;

% get sessions
image_types = fieldnames(lookup_table);

% filenames
for rowIdx = 1:height(table)
    % Allocate based on file name if it's already using BIDS file naming
    if isBIDSFileName([table.Filename{rowIdx}, '.nii'])
        parsed = parseBIDSFilePath([table.Filename{rowIdx}, '.nii']);
        if isfield(parsed, 'ses')
            table_preallocated.Session(rowIdx) = parsed.ses;
        end

        if isfield(parsed, 'suffix')
            table_preallocated.Modality(rowIdx) = parsed.suffix;
            if strcmp(parsed.suffix, 'bold')
                table_preallocated.Type(rowIdx) = 'func';
            elseif strcmp(parsed.suffix, 'dwi')
                table_preallocated.Type(rowIdx) = 'dwi';
            else
                table_preallocated.Type(rowIdx) = 'anat';
            end
        end

        if isfield(parsed, 'run')
            table_preallocated.Run(rowIdx) = parsed.run;
        end

        if isfield(parsed, 'task')
            table_preallocated.Task(rowIdx) = parsed.task;
        end

        if isfield(parsed, 'acq')
            table_preallocated.Acquisition(rowIdx) = parsed.acq;
        end
        continue;
    end

    % Remove folder name from file name when populating the table
    fname = regexprep(table.Filename{rowIdx}, '^sub-[^\W_]+_', '');

    % image types (anat, func, ...)
    for img_type_idx = 1:length(image_types)

        img_type = image_types{img_type_idx};
        sessions = fieldnames(lookup_table.(img_type));

        % sessions (preop, postop, ...)
        for session_idx = 1:length(sessions)

            session = sessions{session_idx};
            modalities = fieldnames(lookup_table.(img_type).(session));

            % modalities (T1w, T2w, ...)
            for modality_idx = 1:length(modalities)

                modality = modalities{modality_idx};
                for name_idx = 1:length(lookup_table.(img_type).(session).(modality))

                    name = char(lookup_table.(image_types{img_type_idx}).(sessions{session_idx}).(modalities{modality_idx}){name_idx});

                    if regexpi(fname, name)
                        table_preallocated.Session(rowIdx) = session;
                        table_preallocated.Type(rowIdx) = img_type;
                        table_preallocated.Modality(rowIdx) = modality;
                    end
                end
            end
        end
    end

    % prepopulate acq tag by resolution for preop MRIs (just anat)
    if  ~strcmp(string(table_preallocated.Session(rowIdx)), 'postop') && ~strcmp(string(table_preallocated.Type(rowIdx)), 'func') && ~strcmp(string(table_preallocated.Type(rowIdx)), 'dwi')
        table_preallocated.Acquisition(rowIdx) = ea_checkacq(imgs_resolution{rowIdx});
    end

end

end

%% cancel button
function cancel_button_function(uiapp)

s = uiconfirm(uiapp.UIFigure, 'Do you really want to cancel file selection?', 'Confirm close', ...
    'Options', {'Yes', 'No'}, 'Icon', 'question');

switch s
    case 'Yes'
        sortedFiles = [];
        setappdata(groot, 'sortedFiles', sortedFiles);
        setappdata(groot, 'returnCode', 'cancel');
        delete(uiapp);
end
end

%% ok button
function ok_button_function(uiapp, table_options, dataset_folder, subjID, postop_modalities, postop_acq_tags)

% sanity checks first
% if preop is empty
nopostop_set = 0;   % we will use this to detect if no postop as been set, because in this case the user may proceed
if isempty(uiapp.previewtree_preop_anat.Children)
    uialert(uiapp.UIFigure, 'No preop files included. Please select at least one preop image.', 'Invalid file selection');
    return

    % if postop is empty, only give out warning
elseif isempty(uiapp.previewtree_postop_anat.Children)
    s = uiconfirm(uiapp.UIFigure, 'No postop anatomical files included. Do you still want to proceed?', 'Confirm missing postop files', ...
        'Options', {'Yes', 'No'}, 'Icon', 'warning');

    switch s
        case 'No'
            return
        case 'Yes'
            nopostop_set = 1;
    end
end

% if multiple files for same session/modality/type are detected
duplicates_found = 0;
N_sessions = length(table_options.Session);
N_types = length(table_options.Type);

for s = 1:N_sessions
    for t = 1:N_types
        uifield = ['previewtree_', table_options.Session{s}, '_', table_options.Type{t}];
        if ~isa(uiapp.(uifield).Children, 'matlab.graphics.GraphicsPlaceholder')
            if contains([uiapp.(uifield).Children.Text], '>>')
                duplicates_found = 1;
            end
        end
    end
end

if duplicates_found
    uialert(uiapp.UIFigure, 'Multiple files with same modality inluded (look for >> filename << in preview window). Please select only one file per modality and session or seperate them by specifying a description.', 'Warning', 'Icon','warning');
    return
end

% go through all the files, check if session, type and modality have been set correctly
postop_modality_found = 0;
postop_mri_found = 0;
for i = find(uiapp.niiFileTable.Data.Include)'
    session = char(uiapp.niiFileTable.Data.Session(i));
    type = char(uiapp.niiFileTable.Data.Type(i));
    modality = char(uiapp.niiFileTable.Data.Modality(i));

    % check whether there is one that has not been defined properly
    if ismember('-', {session, type, modality})
        uialert(uiapp.UIFigure, 'Please specify session, type and modality for all Included images.', 'Invalid file selection');
        return
    end

    % check if postop images have the correct modality
    if strcmp(session, 'postop') && ismember(modality, postop_modalities)
        postop_modality_found = 1;
    end

    if strcmp(session, 'postop') && strcmp(modality, 'MRI')
        postop_mri_found = 1;
    end

end

% if a postop image should be included and has not been set properly, catch this here
if ~(postop_modality_found == 1) && ~(nopostop_set == 1)    % only halt if user has specified postop, but it has the wrong modality
    warning_str = ['No valid modality for the postop session has been found, please choose one of the following:', newline, ...
        sprintf('%s, ', postop_modalities{:})];
    uialert(uiapp.UIFigure, warning_str, 'Invalid file selection');
    return
end

% for the postop MRIs, check whether acquisition tags have been set correctly
if ~(nopostop_set == 1) && postop_mri_found == 1

    ax_postop_mri_found = 0;
    for i = find(uiapp.niiFileTable.Data.Include)'
        session = char(uiapp.niiFileTable.Data.Session(i));
        modality = char(uiapp.niiFileTable.Data.Modality(i));
        acq = char(uiapp.niiFileTable.Data.Acquisition(i));

        % check if postop images have the correct modality
        if strcmp(session, 'postop') && strcmp(modality, 'MRI')
            if ~ismember(acq, postop_acq_tags)
                uialert(uiapp.UIFigure, 'For postop MRIs, the acquisition tag may only be set to <ax>, <sag> or <cor>.', 'Invalid acquisition tag in postop MRI');
                return
            end
            if strcmp(acq, 'ax')
                ax_postop_mri_found = 1;
            end

        end
    end

    if ~(ax_postop_mri_found == 1)
        uialert(uiapp.UIFigure, 'No postop MRI with the acquisition <ax> found. At least one needs to be present that will be used to reconstruct electrode position.', ...
            'Invalid acquisition tag in postop MRI');
        return
    end
end

sortedFiles = cell2struct(cell(1,N_sessions), table_options.Session, N_sessions);

for i = find(uiapp.niiFileTable.Data.Include)'

    session = char(uiapp.niiFileTable.Data.Session(i));
    type = char(uiapp.niiFileTable.Data.Type(i));
    run = char(uiapp.niiFileTable.Data.Run(i));
    task = char(uiapp.niiFileTable.Data.Task(i));
    modality = char(uiapp.niiFileTable.Data.Modality(i));
    acq = char(uiapp.niiFileTable.Data.Acquisition(i));

    % depending on the modality, choose extensions of files to be copied
    if ~strcmp(modality, 'dwi')
        extensions = {'.nii.gz', '.json'};
    else
        extensions = {'.nii.gz', '.json', '.bval', '.bvec'};
    end

    % get filename
    fname = char(uiapp.niiFileTable.UserData.fnames(i));

    export_folder = fullfile(dataset_folder, 'rawdata', subjID, ['ses-', session], type);
    if ~isfolder(export_folder)
        mkdir(export_folder);
    end

    destin_no_ext = fullfile(export_folder, fname);
    source_no_ext = fullfile(getappdata(uiapp.UIFigure, 'niiFolder'), uiapp.niiFileTable.Data.Filename{i});

    for j = 1:length(extensions)
        source = [source_no_ext extensions{j}];
        if isfile(source)
            copyfile(source, [destin_no_ext, extensions{j}])
        else
            ea_cprintf('CmdWinWarnings', 'File not found:\n%s\n', source);
        end
    end

    % generate key for .json file
    if strcmp('-', acq) || strcmp('', acq)
        acq_mod = modality;
    else
        acq_mod = [acq, '_', modality];
    end

    % add file to sortedFiles
    if ~isfield(sortedFiles.(session), type) || ~isfield(sortedFiles.(session).(type), acq_mod) % if no other filename exists for this combination
        sortedFiles.(session).(type).(acq_mod) = {fname}; % set output struct
    else % otherwise, append
        sortedFiles.(session).(type).(acq_mod){end+1} = fname;
    end

end

setappdata(groot, 'sortedFiles', sortedFiles);
setappdata(groot, 'returnCode', 'okay');
delete(uiapp);      % close window

end

%% preview images in the middle
function preview_nii(uiapp, imgs, event)

if isempty(event)
    img_idx = 1;
elseif isempty(event.Indices)
    img_idx = [];
else
    img_idx = event.Indices(1);
end

if ~isempty(img_idx)

    img = imgs{img_idx, 1};

    % update info area
    try
        time_and_date_ugly = img.p.nii.ext.edata_decoded.AcquisitionDateTime;
        time_and_date_pretty = sprintf('%s.%s.%s %s:%s', num2str(time_and_date_ugly(7:8)), ...
            num2str(time_and_date_ugly(5:6)), num2str(time_and_date_ugly(1:4)), ...
            num2str(time_and_date_ugly(9:10)), num2str(time_and_date_ugly(11:12)));
    catch
        time_and_date_pretty = 'N/A';
    end
    info_str = sprintf(['Size:\t\t\t[%s x %s x %s x %s]\n', ...
        'Pixel dimensions:\t[%.2f x %.2f x %.2f]\n', ...
        'Acquistion date:\t%s\n', ...
        'Intensity range:\t[%.0f, %.0f]\n', ...
        'Histogram range:\t[%.0f, %.0f]\n', ...
        'Flip:\t\t\t\t[%d %d %d]'], ...
        num2str(img.p.nii.hdr.dim(2)),num2str(img.p.nii.hdr.dim(3)), num2str(img.p.nii.hdr.dim(4)), num2str(img.p.nii.hdr.dim(5)), ...
        img.p.pixdim(1), img.p.pixdim(2), img.p.pixdim(3), ...
        time_and_date_pretty, ...
        min(img.p.nii.img(:)), max(img.p.nii.img(:)), ...
        min(img.img_thresholded), max(img.img_thresholded), ...
        img.p.flip(1), img.p.flip(2), img.p.flip(3));

    % if .json has been found, insert this into the info string as well
    if img.json_found == 1
        info_str = sprintf('%s\n\nInfo found in JSON sidecar:\n', info_str);

        keys = fieldnames(img.json_sidecar);
        for i = 1:length(keys)
            try
                value = getfield(img.json_sidecar, keys{i});
                if ~ischar(value)
                    value = num2str(value);
                end
                info_str = sprintf('%s\n%s:\t%s', info_str, keys{i}, value);
            end
        end
    end
    uiapp.infoArea.Value = {info_str};

    % update histgram
    h = histogram(uiapp.histogramAxes, img.img_thresholded, 'EdgeAlpha', 0.1, 'FaceColor', [1 1 1], 'EdgeColor', [1 1 1]);
    uiapp.histogramAxes.Color = [0,0,0];

    % plot images
    setappdata(uiapp.UIFigure, 'img', img);

    % coronal
    cut_slice = round(img.dim(2)/2);
    imagesc(uiapp.axes_cor, rot90(squeeze(img.p.nii.img(:, cut_slice, :, 1))), 'ButtonDownFcn', @(src, event) sliceButtonDownFunc(uiapp, event));
    uiapp.axes_cor.Colormap = gray(128);
    setappdata(uiapp.UIFigure, 'cut_slice_cor', cut_slice); % save current cut slice for scrolling
    uiapp.axes_cor.DataAspectRatioMode = 'manual';
    uiapp.axes_cor.DataAspectRatio = [img.p.pixdim(3), img.p.pixdim(1), 1];

    uiapp.axes_cor.YLabel.String = 'L';
    uiapp.axes_cor.YLabel.Color = 'w';
    uiapp.axes_cor.YLabel.Rotation = 0;
    uiapp.axes_cor.YLabel.Position(2) = img.dim(3)/2 + uiapp.axes_cor.YLabel.Extent(4)/2;
    uiapp.axes_cor.YLabel.Position(1) = -3;

    uiapp.axes_cor.Title.String = 'S';
    uiapp.axes_cor.Title.Color = 'w';
    uiapp.axes_cor.Title.Position(1:2) = [img.dim(1)/2, 0];

    % sagittal
    cut_slice = round(img.dim(1)/2);
    imagesc(uiapp.axes_sag, rot90(squeeze(img.p.nii.img(cut_slice, :, :, 1))), 'ButtonDownFcn', @(src, event) sliceButtonDownFunc(uiapp, event));
    uiapp.axes_sag.Colormap = gray(128);
    setappdata(uiapp.UIFigure, 'cut_slice_sag', cut_slice); % save current cut slice for scrolling
    uiapp.axes_sag.DataAspectRatioMode = 'manual';
    uiapp.axes_sag.DataAspectRatio = [img.p.pixdim(3), img.p.pixdim(2), 1];

    uiapp.axes_sag.YLabel.String = 'P';
    uiapp.axes_sag.YLabel.Color = 'w';
    uiapp.axes_sag.YLabel.Rotation = 0;
    uiapp.axes_sag.YLabel.Position(2) = img.dim(3)/2 + uiapp.axes_sag.YLabel.Extent(4)/2;
    uiapp.axes_sag.YLabel.Position(1) = -3;

    uiapp.axes_sag.Title.String = 'S';
    uiapp.axes_sag.Title.Color = 'w';
    uiapp.axes_sag.Title.Position(1:2) = [img.dim(2)/2, 0];

    % axial
    cut_slice = round(img.dim(3)/2);
    imagesc(uiapp.axes_axi, rot90(img.p.nii.img(:, :, cut_slice)), 'ButtonDownFcn', @(src, event) sliceButtonDownFunc(uiapp, event));
    uiapp.axes_axi.Colormap = gray(128);
    setappdata(uiapp.UIFigure, 'cut_slice_axi', cut_slice); % save current cut slice for scrolling
    uiapp.axes_axi.DataAspectRatioMode = 'manual';
    uiapp.axes_axi.DataAspectRatio = [img.p.pixdim(2), img.p.pixdim(1), 1];

    uiapp.axes_axi.YLabel.String = 'L';
    uiapp.axes_axi.YLabel.Color = 'w';
    uiapp.axes_axi.YLabel.Rotation = 0;
    uiapp.axes_axi.YLabel.Position(2) = img.dim(2)/2 + uiapp.axes_axi.YLabel.Extent(4)/2;
    uiapp.axes_axi.YLabel.Position(1) = -3;

    uiapp.axes_axi.Title.String = 'S';
    uiapp.axes_axi.Title.Color = 'w';
    uiapp.axes_axi.Title.Position(1:2) = [img.dim(1)/2, 0];

    update_crosschairs(uiapp, img.dim);
end
end

%% scroll images
function scroll_nii(uiapp, event)

axesTag = getMousePointerAxes(uiapp.UIFigure, uiapp.RightPanel);
img = getappdata(uiapp.UIFigure, 'img');
dim = img.dim;

if ~isempty(axesTag)
    sliceUpdated.cor = 0;
    sliceUpdated.sag = 0;
    sliceUpdated.axi = 0;
    switch axesTag
        case 'axi'
            sliceNr = getappdata(uiapp.UIFigure, 'cut_slice_axi');
            if event.VerticalScrollCount == -1 % up scroll
                if ~(sliceNr >= img.dim(3) - 2)
                    sliceNr = sliceNr + 2;
                end
            else
                if ~(sliceNr <= 2)
                    sliceNr = sliceNr - 2;
                end

            end
            imagesc(uiapp.axes_axi, rot90(img.p.nii.img(:, :, sliceNr)), 'ButtonDownFcn', @(src, event) sliceButtonDownFunc(uiapp, event));
            setappdata(uiapp.UIFigure, 'cut_slice_axi', sliceNr);
            sliceUpdated.axi = 1;
        case 'cor'
            sliceNr = getappdata(uiapp.UIFigure, 'cut_slice_cor');
            if event.VerticalScrollCount == -1 % up scroll
                if ~(sliceNr >= img.dim(2) - 2)
                    sliceNr = sliceNr + 2;
                end
            else
                if ~(sliceNr <= 2)
                    sliceNr = sliceNr - 2;
                end

            end
            imagesc(uiapp.axes_cor, rot90(squeeze(img.p.nii.img(:, sliceNr, :, 1))), 'ButtonDownFcn', @(src, event) sliceButtonDownFunc(uiapp, event));
            setappdata(uiapp.UIFigure, 'cut_slice_cor', sliceNr);
            sliceUpdated.cor = 1;
        case 'sag'
            sliceNr = getappdata(uiapp.UIFigure, 'cut_slice_sag');
            if event.VerticalScrollCount == -1 % up scroll
                if ~(sliceNr >= img.dim(1) - 2)
                    sliceNr = sliceNr + 2;
                end
            else
                if ~(sliceNr <= 2)
                    sliceNr = sliceNr - 2;
                end

            end
            imagesc(uiapp.axes_sag, rot90(squeeze(img.p.nii.img(sliceNr, :, :, 1))), 'ButtonDownFcn', @(src, event) sliceButtonDownFunc(uiapp, event));
            setappdata(uiapp.UIFigure, 'cut_slice_sag', sliceNr);
            sliceUpdated.sag = 1;
    end
    update_crosschairs(uiapp, dim, sliceUpdated);
end

end

%% update cross
function update_crosschairs(uiapp, dim, sliceUpdated)
    corSliceNr = getappdata(uiapp.UIFigure, 'cut_slice_cor'); % y, dim(2)
    sagSliceNr = getappdata(uiapp.UIFigure, 'cut_slice_sag'); % x, dim(1)
    axiSliceNr = getappdata(uiapp.UIFigure, 'cut_slice_axi'); % z, dim(3)
    
    if ~exist('sliceUpdated', 'var')
        sliceUpdated.cor = 0;
        sliceUpdated.sag = 0;
        sliceUpdated.axi = 0;
    end

    corAxesHorzLine = findobj(uiapp.axes_cor, 'Type', 'Line', 'Tag', 'HorzLine'); % horz line is axi slice
    if isempty(corAxesHorzLine) || sliceUpdated.axi
        if sliceUpdated.axi
            delete(corAxesHorzLine);
        end
        line(uiapp.axes_cor, [1, dim(1)], [dim(3)-axiSliceNr, dim(3)-axiSliceNr], 'Color', 'b', 'LineWidth', 1, 'Tag', 'HorzLine');
    end

    corAxesVertLine = findobj(uiapp.axes_cor, 'Type', 'Line', 'Tag', 'VertLine'); % vert line is sag slice
    if isempty(corAxesVertLine) || sliceUpdated.sag
        if sliceUpdated.sag
            delete(corAxesVertLine);
        end
        line(uiapp.axes_cor, [sagSliceNr, sagSliceNr], [1, dim(3)], 'Color', 'b', 'LineWidth', 1, 'Tag', 'VertLine');
    end

    sagAxesHorzLine = findobj(uiapp.axes_sag, 'Type', 'Line', 'Tag', 'HorzLine'); % horz line is axi slice
    if isempty(sagAxesHorzLine) || sliceUpdated.axi
        if sliceUpdated.axi
            delete(sagAxesHorzLine);
        end
        line(uiapp.axes_sag, [1, dim(2)], [dim(3)-axiSliceNr, dim(3)-axiSliceNr], 'Color', 'b', 'LineWidth', 1, 'Tag', 'HorzLine');
    end

    sagAxesVertLine = findobj(uiapp.axes_sag, 'Type', 'Line', 'Tag', 'VertLine'); % vert line is cor slice
    if isempty(sagAxesVertLine) || sliceUpdated.cor
        if sliceUpdated.cor
            delete(sagAxesVertLine);
        end
        line(uiapp.axes_sag, [corSliceNr, corSliceNr], [1, dim(3)], 'Color', 'b', 'LineWidth', 1, 'Tag', 'VertLine');
    end

    axiAxesHorzLine = findobj(uiapp.axes_axi, 'Type', 'Line', 'Tag', 'HorzLine'); % horz line is cor slice
    if isempty(axiAxesHorzLine) || sliceUpdated.cor
        if sliceUpdated.cor
            delete(axiAxesHorzLine);
        end
        line(uiapp.axes_axi, [1, dim(1)], [dim(2)-corSliceNr, dim(2)-corSliceNr], 'Color', 'b', 'LineWidth', 1, 'Tag', 'HorzLine');
    end

    axiAxesVertLine = findobj(uiapp.axes_axi, 'Type', 'Line', 'Tag', 'VertLine'); % vert line is sag slice
    if isempty(axiAxesVertLine) || sliceUpdated.sag
        if sliceUpdated.sag
            delete(axiAxesVertLine);
        end
        line(uiapp.axes_axi, [sagSliceNr, sagSliceNr], [1, dim(2)], 'Color', 'b', 'LineWidth', 1, 'Tag', 'VertLine');
    end
end

%% Button down on axes
function sliceButtonDownFunc(uiapp, event)
    img = getappdata(uiapp.UIFigure, 'img');
    dim = img.dim;
    if ~isempty(img)
        x = round(event.IntersectionPoint(1));
        y = round(event.IntersectionPoint(2));
        sliceUpdated.cor = 0;
        sliceUpdated.sag = 0;
        sliceUpdated.axi = 0;
        switch event.Source.Parent.Tag
            case 'cor'
                imagesc(uiapp.axes_sag, rot90(squeeze(img.p.nii.img(x, :, :))), 'ButtonDownFcn', @(src, event) sliceButtonDownFunc(uiapp, event));
                setappdata(uiapp.UIFigure, 'cut_slice_sag', x);
                sliceUpdated.sag = 1;
                imagesc(uiapp.axes_axi, rot90(img.p.nii.img(:, :, dim(3)-y)), 'ButtonDownFcn', @(src, event) sliceButtonDownFunc(uiapp, event));
                setappdata(uiapp.UIFigure, 'cut_slice_axi', dim(3)-y);
                sliceUpdated.axi = 1;
            case 'sag'
                imagesc(uiapp.axes_cor, rot90(squeeze(img.p.nii.img(:, x, :))), 'ButtonDownFcn', @(src, event) sliceButtonDownFunc(uiapp, event));
                setappdata(uiapp.UIFigure, 'cut_slice_cor', x);
                sliceUpdated.cor = 1;
                imagesc(uiapp.axes_axi, rot90(img.p.nii.img(:, :, dim(3)-y)), 'ButtonDownFcn', @(src, event) sliceButtonDownFunc(uiapp, event));
                setappdata(uiapp.UIFigure, 'cut_slice_axi', dim(3)-y);
                sliceUpdated.axi = 1;
            case 'axi'
                imagesc(uiapp.axes_sag, rot90(squeeze(img.p.nii.img(x, :, :))), 'ButtonDownFcn', @(src, event) sliceButtonDownFunc(uiapp, event));
                setappdata(uiapp.UIFigure, 'cut_slice_sag', x);
                sliceUpdated.sag = 1;
                imagesc(uiapp.axes_cor, rot90(squeeze(img.p.nii.img(:, dim(2)-y, :))), 'ButtonDownFcn', @(src, event) sliceButtonDownFunc(uiapp, event));
                setappdata(uiapp.UIFigure, 'cut_slice_cor', dim(2)-y);
                sliceUpdated.cor = 1;
        end
        update_crosschairs(uiapp, dim, sliceUpdated);
    end
end

%% check where the mouse pointer is
function axesTag = getMousePointerAxes(fig, panel)

oldUnits = get(0,'units');
set(0,'units','pixels');
% Get the figure beneath the mouse pointer & mouse pointer pos

p = get(0,'PointerLocation');
set(0,'units',oldUnits);
% Look for quick exit (if mouse pointer is not over any figure)

% Compute figure offset of mouse pointer in pixels
figPos = getpixelposition(fig);
panelPos = getpixelposition(panel);

% Loop over all figure descendants
c = findobj(get(fig,'Children'), 'type', 'axes');
for h = c'

    axesPos = getpixelposition(h);  % Note: cache this for improved performance

    x_lower = figPos(1) + panelPos(1) + axesPos(1);
    x_upper = x_lower + axesPos(3);

    y_lower = figPos(2) + panelPos(2) + axesPos(2);
    y_upper = y_lower + axesPos(4);

    % If descendant contains the mouse pointer position, exit

    if (p(1) > x_lower) && (p(1) < x_upper) && (p(2) > y_lower) && (p(2) < y_upper)
        axesTag = h.Tag;
        return
    end
end
axesTag = [];
end


%% helper functions for lookup_table
function lookup_table_data = convert_lookup_struct_to_table(lookup_table)
types = fieldnames(lookup_table);

lookup_table_data = table();

for type_idx = 1:length(types)

    type = types{type_idx};

    sessions = fieldnames(lookup_table.(type));

    for ses_idx = 1:length(sessions)

        session = sessions{ses_idx};

        modalities = fieldnames(lookup_table.(type).(session));

        temp_table = table();
        for mod_idx = 1:length(modalities)

            modality = modalities{mod_idx};

            keywords =  lookup_table.(type).(session).(modality);
            keywords_table = "";

            if ~isempty(keywords)

                for keyword_idx = 1:length(keywords)
                    if keyword_idx == 1
                        keywords_table = string(sprintf('%s', keywords{keyword_idx}));
                    else
                        keywords_table = string(sprintf('%s,%s', keywords_table, keywords{keyword_idx}));
                    end
                end
            end

            temp_table.Keywords = string(keywords_table);
            temp_table.Type = string(type);
            temp_table.Session = string(session);
            temp_table.Modality = string(modality);

            lookup_table_data = [lookup_table_data; temp_table];
        end
    end
end

end


%% misc helper functions
function bids_fname = generate_bids_filename(subjID, session, run, task, acq, modality, phase_enc_dir)

tag_names = {{'ses' session}, {'run', run}, {'task', task}, {'acq', acq}};

bids_fname = sprintf('%s', subjID);

for i = 1:numel(tag_names)
    if ~(strcmp(tag_names{1, i}{2}, '-')) && ~(strcmp(tag_names{1, i}{2}, '')) && ~isempty(tag_names{1, i}{2})
        bids_fname = sprintf('%s_%s-%s', bids_fname, tag_names{1, i}{1}, tag_names{1, i}{2});
    end
end

if ~isempty(phase_enc_dir)
    bids_fname = [bids_fname, '_', 'dir-', phase_enc_dir, '_', modality];
else
    bids_fname = [bids_fname, '_', modality];
end

end


function ped = get_phase_encoding_direction(json_sidecar)

if isfield(json_sidecar, 'PhaseEncodingDirection')
    switch json_sidecar.PhaseEncodingDirection
        case 'j-'
            ped = 'AP';
        case 'j'
            ped = 'PA';
        case 'i'
            ped = 'LR';
        case 'i-'
            ped = 'RL';
        otherwise
            ped = '';
    end
else
    ped = '';
end
end


function dupl_row = find_dupl_file(niitable, current_row, session, type, run, task, modality)

dupl_row = [];
for i = 1:height(niitable.Data)
    if strcmp(session, char(niitable.Data.Session(i))) ...
            && strcmp(type, char(niitable.Data.Type(i))) ...
            && strcmp(run, char(niitable.Data.Run(i))) ...
            && strcmp(task, char(niitable.Data.Task(i))) ...
            && strcmp(modality, char(niitable.Data.Modality(i))) ...
            && i ~= current_row ...
            && niitable.Data.Include(i) == 1
        dupl_row = i;
    end
end
end
