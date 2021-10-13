function anat_files = ea_dicom_to_bids(subjID, fnames, dataset_folder)

% Function creates a GUI in order to select which files should be used to create a BIDS compliant dataset and
% which files should be used by lead-dbs. After selection, a rawdata folder will be created and selected files will be
% copy/pasted to the appropriate locations
%
% input
%   subjID (string): subjID, including the 'sub' tag
%   fnames (cell): cell of filenames to be loaded
%   dataset_folder (string): root folder of BIDS-like dataset
%
% The script expects a sourcedata/subjID/ folder with DICOM files to be present under the dataset_folder path
%
% output
%   anat_files (struct): struct with preop and postop session fields, inside each session field there are fields for every modality
%                                 those fields have strings with the filenames that are selected and to be used by lead-dbs
% __________________________________________________________________________________
% Copyright (C) 2021 Charite University Medicine Berlin, Movement Disorders Unit
% Johannes Achtzehn
% Part of this code is inspired by the dicm2nii tool of xiangrui.li@gmail.com

% lookup table to pre-allocate options for the table.% check whether pref exists, if not create it
if ~ispref('dcm2bids', 'lookuptable')
    lookup_table = loadjson(fullfile(ea_getearoot(), 'helpers', 'dicom_bids_lookuptable.json'));
    addpref('dcm2bids', 'lookuptable', lookup_table)
else
    lookup_table = getpref('dcm2bids', 'lookuptable');
end

nii_folder = fullfile(dataset_folder, 'sourcedata', subjID, 'tmp');    % where are the nifti files located?
N_fnames = length(fnames);

imgs = cell(N_fnames,1);
h_wait = waitbar(0, 'Please wait while Niftii images are being loaded');
for image_idx = 1:N_fnames
    imgs{image_idx} = struct();
    [imgs{image_idx}.p, imgs{image_idx}.frm, imgs{image_idx}.rg, imgs{image_idx}.dim] = read_nii(fullfile(nii_folder, [fnames{image_idx}, '.nii.gz']), [], 0);
    imgs{image_idx}.percentile = prctile(imgs{image_idx}.p.nii.img(:), 95, 'all');
    imgs{image_idx}.img_thresholded = imgs{image_idx}.p.nii.img(:);
    imgs{image_idx}.img_thresholded(imgs{image_idx}.img_thresholded > imgs{image_idx}.percentile) = nan;
    imgs{image_idx}.img_thresholded(imgs{image_idx}.img_thresholded < median(imgs{image_idx}.p.nii.img(:))) = nan;
    imgs{image_idx}.img_thresholded = imgs{image_idx}.img_thresholded(~isnan(imgs{image_idx}.img_thresholded));
    
    waitbar(image_idx / N_fnames, h_wait, sprintf('Please wait while Niftii images are being loaded (%i/%i)', image_idx, length(fnames)));
end
close(h_wait);

% options that should appear in the table
table_options = struct;
table_options.Session = {'preop', 'postop'};
table_options.Type = {'anat'};
table_options.Modality = {'T1w', 'T2w', 'FGATIR', 'FLAIR', 'T2starw', 'PDw', 'CT', 'ax_MR', 'sag_MR', 'cor_MR'};

nifti_table = structfun(@(x) categorical(repmat({'-'}, [N_fnames,1]), ['-' x]), table_options, 'uni', 0);
nifti_table.Include = false(N_fnames,1);
nifti_table.Filename = fnames;

nifti_table = orderfields(nifti_table, [4 5 1 2 3]);
nifti_table = struct2table(nifti_table);

try
    nifti_table_preallocated = preallocate_table(nifti_table, lookup_table);
catch
    disp('Preallocation of table failed, defaulting to skip!');
    nifti_table_preallocated = nifti_table;
end

% create GUI
uiapp = dicom_to_bids;

% populate table
uiapp.niiFileTable.Data = nifti_table_preallocated;
uiapp.niiFileTable.ColumnEditable = [true false true true true];

% set subject ID and file path
uiapp.SubjectIDLabel.Text = sprintf('Conversion for subject: %s',subjID);
uiapp.FilepathLabel.Text = sprintf('BIDS root: %s',dataset_folder);

% update preview tree and expand it
uiapp.previewtree_subj.Text = subjID;
expand(uiapp.Tree, 'all');

preview_nii(uiapp, imgs{1,1}); % set initial image to the first one
cell_change_callback(uiapp, table_options, subjID, []) % call preview tree updater to get preallocated changes

uiapp.niiFileTable.CellSelectionCallback = @(src,event) preview_nii(uiapp,imgs{event.Indices(1), 1}); % callback for table selection -> display current selected image
uiapp.niiFileTable.CellEditCallback = @(src,event) cell_change_callback(uiapp, table_options, subjID, event); % callback for cell change -> update uiapp tree on the right

uiapp.UIFigure.WindowScrollWheelFcn = @(src, event) scroll_nii(uiapp, event);     % callback for scrolling images

% OK button behaviour
uiapp.OKButton.ButtonPushedFcn = @(btn,event) ok_button_function(uiapp, table_options, dataset_folder, nii_folder, subjID);

% cancel button behaviour
uiapp.CancelButton.ButtonPushedFcn =  @(btn,event) cancel_button_function(uiapp);

% looup table behaviour
uiapp.LookupButton.ButtonPushedFcn = @(btn,event) lookup_button_function(uiapp);
waitfor(uiapp.UIFigure);

try
    anat_files = getappdata(groot, 'anat_files');
catch
    anat_files = [];
    
end

end

%% lookup button function
function lookup_button_function(uiapp)

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
lookup_table_gui.SaveButton.ButtonPushedFcn = @(btn,event) save_lookup_function(uiapp, lookup_table_gui);

waitfor(lookup_table_gui.UIFigure);
end

function save_lookup_function(main_gui, lookup_table_gui)

lookup_table = convert_table_to_lookup_struct(lookup_table_gui.UITable.Data);

setpref('dcm2bids', 'lookuptable', lookup_table);

T_preallocated = preallocate_table(main_gui.niiFileTable.Data, lookup_table);

main_gui.niiFileTable.Data = T_preallocated;

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
function cell_change_callback(uiapp, table_options, subjID, event)

uiapp.previewtree_preop_anat.Children.delete;      % delete children
uiapp.previewtree_postop_anat.Children.delete;    % delete children

% go through all the files that have been selected to include
for i = find(uiapp.niiFileTable.Data.Include)'
    
    session = char(uiapp.niiFileTable.Data.Session(i));
    type = char(uiapp.niiFileTable.Data.Type(i));
    modality = char(uiapp.niiFileTable.Data.Modality(i));
 
    % check whether everything has been properly defined befor updating uitree
    if ~any(strcmp('-', {session, type, modality}))
        fname = sprintf('%s_ses-%s_%s', subjID, session, modality);   % generate BIDS filename
        ui_field = ['previewtree_' session '_anat'];
        if ~isempty(uiapp.(ui_field).Children) && any(ismember(fname, {uiapp.(ui_field).Children.Text}))
            fname = ['>> ', fname, ' <<'];
        end
        uitreenode(uiapp.(ui_field), 'Text', fname);
    end
    
end

for i = find(~uiapp.niiFileTable.Data.Include)'
     session = char(uiapp.niiFileTable.Data.Session(i));
    type = char(uiapp.niiFileTable.Data.Type(i));
    modality = char(uiapp.niiFileTable.Data.Modality(i));
    
    if ~isempty(event)
    if ~any(strcmp('-', {session, type, modality})) && event.Indices(2) > 2
        uiapp.niiFileTable.Data.Include(i) = true;
        fname = sprintf('%s_ses-%s_%s', subjID, session, modality);   % generate BIDS filename
        ui_field = ['previewtree_' session '_anat'];
        if ~isempty(uiapp.(ui_field).Children) && any(ismember(fname, {uiapp.(ui_field).Children.Text}))
            fname = ['>> ', fname, ' <<'];
        end
        uitreenode(uiapp.(ui_field), 'Text', fname);
    end
    end
    
end
    
end

%% preallocate table on the left
function table_preallocated = preallocate_table(table, lookup_table)

table_preallocated = table;

% get sessions
image_types = fieldnames(lookup_table);

% filenames
for rowIdx = 1:height(table)
    
    fname = table.Filename{rowIdx};
    
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
                    
                    if regexp(fname, name)
                        table_preallocated.Session(rowIdx) = session;
                        table_preallocated.Type(rowIdx) = img_type;
                        table_preallocated.Modality(rowIdx) = modality;
                        table_preallocated.Include(rowIdx) = true;
                    end
                end
            end
        end
    end
end

end

%% cancel button
function cancel_button_function(uiapp)

s = uiconfirm(uiapp.UIFigure, 'Do you really want to cancel file selection?', 'Confirm close', ...
    'Options', {'Yes', 'No'}, 'Icon', 'question');

switch s
    case 'Yes'
        delete(uiapp);
end


end

%% ok button
function ok_button_function(uiapp, table_options, dataset_folder, nii_folder, subjID)

% sanity check
if isempty(uiapp.previewtree_preop_anat.Children)
    uialert(uiapp.UIFigure, 'No preop files included. Please select at least one preop image.', 'Warning', 'Icon','warning');
    return
elseif contains([uiapp.previewtree_preop_anat.Children.Text], '>>') || ...
        (~isempty(uiapp.previewtree_postop_anat.Children) && contains([uiapp.previewtree_postop_anat.Children.Text], '>>'))
    uialert(uiapp.UIFigure, 'Multiple files with same modality inluded (look for >> filename << in preview window). Please select only one file per modality and session.', 'Warning', 'Icon','warning');
    return
else
    for i = find(uiapp.niiFileTable.Data.Include)'
        session = char(uiapp.niiFileTable.Data.Session(i));
    type = char(uiapp.niiFileTable.Data.Type(i));
    modality = char(uiapp.niiFileTable.Data.Modality(i));
    include = uiapp.niiFileTable.Data.Include(i);
    
    % check whether there is one that has not been defined properly
    if any(strcmp('-', {session, type, modality}))
        uialert(uiapp.UIFigure, 'Please specify session, type and modality for all Included images.', 'Warning', 'Icon','warning');
       return
    end
    end
    
end

N_sessions = length(table_options.Session);
anat_files = cell2struct(cell(1,N_sessions), table_options.Session, N_sessions);

extensions = {'.nii.gz', '.json'};

for i = find(uiapp.niiFileTable.Data.Include)'
    
    session = char(uiapp.niiFileTable.Data.Session(i));
    modality = char(uiapp.niiFileTable.Data.Modality(i));
    
    fname = sprintf('%s_ses-%s_%s', subjID, session, modality);
    
    export_folder = fullfile(dataset_folder, 'rawdata', subjID, ['ses-', session], 'anat');
    if ~isfolder(export_folder)
        mkdir(export_folder);
    end
    
    destin_no_ext = fullfile(export_folder, fname);
    source_no_ext = fullfile(nii_folder, uiapp.niiFileTable.Data.Filename{i});
    
    for j = 1:length(extensions)
        source = [source_no_ext extensions{j}];
        if isfile(source)
            copyfile(source, [destin_no_ext extensions{j}])
        else
            warning('Selected file %s cannot be found and was not copied! Please copy/paste manually.\n', source);
        end
    end
    
    anat_files.(session).(modality) = fname; % set output struct
    
end

setappdata(groot, 'anat_files', anat_files);
delete(uiapp);      % close window

end

%% preview images in the middle
function preview_nii(uiapp, img)

% update info area
try
    time_and_date_ugly = img.p.nii.ext.edata_decoded.AcquisitionDateTime;
    time_and_date_pretty = sprintf('%s.%s.%s %s:%s', num2str(time_and_date_ugly(7:8)), ...
        num2str(time_and_date_ugly(5:6)), num2str(time_and_date_ugly(1:4)), ...
        num2str(time_and_date_ugly(9:10)), num2str(time_and_date_ugly(11:12)));
catch
    time_and_date_pretty = 'N/A';
end
info_str = sprintf('Size:\t\t\t[%s x %s x %s]\nPixel dimensions:\t[%.2f x %.2f x %.2f]\nAcquistion date:\t%s\nIntensity range:\t[%.0f, %.0f]\nHistogram range:\t[%.0f, %.0f]', ...
    num2str(img.dim(1)), num2str(img.dim(2)), num2str(img.dim(3)), ...
    img.p.pixdim(1), img.p.pixdim(2), img.p.pixdim(3), ...
    time_and_date_pretty, ...
    min(img.p.nii.img(:)), max(img.p.nii.img(:)), ...
    min(img.img_thresholded), max(img.img_thresholded));

uiapp.infoArea.Value = {info_str};

% update histgram
h = histogram(uiapp.histogramAxes, img.img_thresholded, 'EdgeAlpha', 0.1, 'FaceColor', [1 1 1], 'EdgeColor', [1 1 1]);
uiapp.histogramAxes.Color = [0,0,0];

% plot images
setappdata(uiapp.UIFigure, 'img', img);

% axial
cut_slice = round(img.dim(3)/2);
imagesc(uiapp.axes_axi, img.p.nii.img(:, :, cut_slice));
uiapp.axes_axi.Colormap = gray(128);
setappdata(uiapp.UIFigure, 'cut_slice_axi', cut_slice); % save current cut slice for scrolling
uiapp.axes_axi.DataAspectRatioMode = 'manual';
uiapp.axes_axi.DataAspectRatio = [img.p.pixdim(1), img.p.pixdim(2), 1];
set(uiapp.axes_axi, 'view', [90, -90]);

% coronal
cut_slice = round(img.dim(2)/2);
imagesc(uiapp.axes_cor, squeeze(img.p.nii.img(:, cut_slice, :)));
uiapp.axes_cor.Colormap = gray(128);
setappdata(uiapp.UIFigure, 'cut_slice_cor', cut_slice); % save current cut slice for scrolling
uiapp.axes_cor.DataAspectRatioMode = 'manual';
uiapp.axes_cor.DataAspectRatio = [img.p.pixdim(1), img.p.pixdim(3), 1];
set(uiapp.axes_cor, 'view', [90, -90]);

% sagittal
cut_slice = round(img.dim(1)/2);
imagesc(uiapp.axes_sag, squeeze(img.p.nii.img(cut_slice, :, :)));
uiapp.axes_sag.Colormap = gray(128);
setappdata(uiapp.UIFigure, 'cut_slice_sag', cut_slice); % save current cut slice for scrolling
uiapp.axes_sag.DataAspectRatioMode = 'manual';
uiapp.axes_sag.DataAspectRatio = [img.p.pixdim(1), img.p.pixdim(3), 1];
set(uiapp.axes_sag, 'view', [90, -90]);

end

%% scroll images
function scroll_nii(uiapp, event)

hAxes = checkMousePointer(uiapp.UIFigure, uiapp.RightPanel);
img = getappdata(uiapp.UIFigure, 'img');
dim = img.dim;

if ~isempty(hAxes)
    switch hAxes.Tag
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
            imagesc(uiapp.axes_axi, img.p.nii.img(:, :, sliceNr));
            setappdata(uiapp.UIFigure, 'cut_slice_axi', sliceNr);
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
            imagesc(uiapp.axes_cor, squeeze(img.p.nii.img(:, sliceNr, :)));
            setappdata(uiapp.UIFigure, 'cut_slice_cor', sliceNr);
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
            imagesc(uiapp.axes_sag, squeeze(img.p.nii.img(sliceNr, :, :)));
            setappdata(uiapp.UIFigure, 'cut_slice_sag', sliceNr);
        otherwise
            
    end
end

end

%% check where the mouse pointer is
function h = checkMousePointer(fig, panel)

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
        return
    end
end
h = [];
end

function [p, frm, rg, dim] = read_nii(fname, ask_code, reOri)
if nargin<2, ask_code = []; end
if ischar(fname), p.nii = nii_tool('load', fname);
else, p.nii = fname; fname = p.nii.hdr.file_name;
end
p.hdr0 = p.nii.hdr; % original hdr
c = p.nii.hdr.intent_code;
if c>=3000 && c<=3099 && isfield(p.nii, 'ext') && any([p.nii.ext.ecode] == 32)
    p.nii = cii2nii(p.nii);
end

if nargin<3 || reOri
    [p.nii, p.perm, p.flip] = nii_reorient(p.nii, 0, ask_code);
else
    p.perm = 1:3;
    p.flip = false(1,3);
end

dim = p.nii.hdr.dim(2:8);
dim(dim<1 | mod(dim,1)~=0) = 1;
if p.nii.hdr.dim(1)>4 % 4+ dim, put all into dim4
    if sum(dim(4:7)>1)>1
        warndlg([fname ' has 5 or more dimension. Dimension above 4 are ' ...
            'all treated as volumes for visualization']);
    end
    dim(4) = prod(dim(4:7)); dim(5:7) = 1;
    p.nii.img = reshape(p.nii.img, [dim size(p.nii.img, 8)]);
end

[p.R, frm] = nii_xform_mat(p.nii.hdr, ask_code);
dim = dim(1:3);
p.pixdim = p.nii.hdr.pixdim(2:4);

if size(p.nii.img,4)<4 && ~isfloat(p.nii.img)
    p.nii.img = single(p.nii.img);
end
if p.nii.hdr.scl_slope==0, p.nii.hdr.scl_slope = 1; end
if p.nii.hdr.scl_slope~=1 || p.nii.hdr.scl_inter~=0
    if isfloat(p.nii.img)
        p.nii.img = p.nii.img * p.nii.hdr.scl_slope + p.nii.hdr.scl_inter;
    else
        p.scl_slope = p.nii.hdr.scl_slope;
        p.scl_inter = p.nii.hdr.scl_inter;
    end
end

% check if ROI labels available: the same file name with .txt extension
if c == 1002 % Label
    [pth, nam, ext] = fileparts(p.nii.hdr.file_name);
    nam1 = fullfile(pth, [nam '.txt']);
    if strcmpi(ext, '.gz') && ~exist(nam1, 'file')
        [~, nam] = fileparts(nam);
        nam1 = fullfile(pth, [nam '.txt']);
    end
    if exist(nam1, 'file') % each line format: 1 ROI_1
        fid = fopen(nam1);
        while 1
            ln = fgetl(fid);
            if ~ischar(ln), break; end
            [ind, a] = strtok(ln);
            ind = str2double(ind);
            try p.labels{ind} = strtrim(a); catch, end
        end
        fclose(fid);
    end
end
rg = get_range(p.nii, isfield(p, 'labels'));
try, p.map = p.nii.NamedMap{1}.map; end
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

%% helper functions for nii reading
function [nii, perm, flp] = nii_reorient(nii, leftHand, ask_code)
if nargin<3, ask_code = []; end
[R, frm] = nii_xform_mat(nii.hdr, ask_code);
dim = nii.hdr.dim(2:4);
pixdim = nii.hdr.pixdim(2:4);
[R, perm, flp] = reorient(R, dim, leftHand);
fps = bitand(nii.hdr.dim_info, [3 12 48]) ./ [1 4 16];
if ~isequal(perm, 1:3)
    nii.hdr.dim(2:4) = dim(perm);
    nii.hdr.pixdim(2:4) = pixdim(perm);
    nii.hdr.dim_info = [1 4 16] * fps(perm)' + bitand(nii.hdr.dim_info, 192);
    nii.img = permute(nii.img, [perm 4:8]);
end
sc = nii.hdr.slice_code;
if sc>0 && flp(fps==3)
    nii.hdr.slice_code = sc+mod(sc,2)*2-1; % 1<->2, 3<->4, 5<->6
end
if isequal(perm, 1:3) && ~any(flp), return; end
if frm(1) == nii.hdr.sform_code % only update matching form
    nii.hdr.srow_x = R(1,:);
    nii.hdr.srow_y = R(2,:);
    nii.hdr.srow_z = R(3,:);
end
if frm(1) == nii.hdr.qform_code
    nii.hdr.qoffset_x = R(1,4);
    nii.hdr.qoffset_y = R(2,4);
    nii.hdr.qoffset_z = R(3,4);
    R0 = normc(R(1:3, 1:3));
    dcm2quat = dicm2nii('', 'dcm2quat', 'func_handle');
    [q, nii.hdr.pixdim(1)] = dcm2quat(R0);
    nii.hdr.quatern_b = q(2);
    nii.hdr.quatern_c = q(3);
    nii.hdr.quatern_d = q(4);
end
for i = find(flp), nii.img = flip(nii.img, i); end


end

function [R, frm] = nii_xform_mat(hdr, ask_code)
% [R, form] = nii_xform_mat(hdr, asked_code);
% Return the transformation matrix from a NIfTI hdr. By default, this returns
% the sform if available. If the optional second input, required form code, is
% provided, this will try to return matrix for that form code. The second
% optional output is the form code of the actually returned matrix.
fs = [hdr.sform_code hdr.qform_code]; % sform preferred
if fs(1)==fs(2), fs = fs(1); end % sform if both are the same
f = fs(fs>=1 & fs<=4); % 1/2/3/4 only
if isempty(f) || ~strncmp(hdr.magic, 'n', 1) % treat it as Analyze
    frm = 0;
    try % try spm style Analyze
        [pth, nam, ext] = fileparts(hdr.file_name);
        if strcmpi(ext, '.gz'), [~, nam] = fileparts(nam); end
        R = load(fullfile(pth, [nam '.mat']));
        R = R.M;
    catch % make up R for Analyze: suppose xyz order with left storage
        R = [diag(hdr.pixdim(2:4)) -(hdr.dim(2:4).*hdr.pixdim(2:4)/2)'; 0 0 0 1];
        R(1,:) = -R(1,:); % use left handed
    end
    return;
end

if numel(f)==1 || nargin<2 || isempty(ask_code) % only 1 avail or no ask_code
    frm = f;
else % numel(f) is 2, numel(ask_code) can be 1 or 2
    frm = f(f == ask_code(1));
    if isempty(frm) && numel(ask_code)>1, frm = f(f == ask_code(2)); end
    if isempty(frm) && any(f==2), frm = 2; end % use confusing code 2
    if isempty(frm), frm = f(1); end % no match to ask_code, use sform
end

if frm(1) == fs(1) % match sform_code or no match
    R = [hdr.srow_x; hdr.srow_y; hdr.srow_z; 0 0 0 1];
else % match qform_code
    R = quat2R(hdr);
end
end

function [R, perm, flp] = reorient(R, dim, leftHand)
% [R, perm, flip] = reorient(R, dim, leftHand)
% Re-orient transformation matrix R (4x4), so it will be diagonal major and
% positive at diagonal, unless the optional third input is true, which requires
% left-handed matrix, where R(1,1) will be negative.
% The second input is the img space dimension (1x3).
% The perm output, like [1 2 3] or a permutation of it, indicates if input R was
% permuted for 3 axis. The third output, flip (1x3 logical), indicates an axis
% (AFTER perm) is flipped if true.
a = abs(R(1:3,1:3));
[~, ixyz] = max(a);
if ixyz(2) == ixyz(1), a(ixyz(2),2) = 0; [~, ixyz(2)] = max(a(:,2)); end
if any(ixyz(3) == ixyz(1:2)), ixyz(3) = setdiff(1:3, ixyz(1:2)); end
[~, perm] = sort(ixyz);
R(:,1:3) = R(:,perm);
flp = R([1 6 11]) < 0; % diag(R(1:3, 1:3))
if nargin>2 && leftHand, flp(1) = ~flp(1); end
rotM = diag([1-flp*2 1]);
rotM(1:3, 4) = (dim(perm)-1) .* flp; % 0 or dim-1
R = R / rotM; % xform matrix after flip
end

%  Estimate lower and upper bound of img display
function rg = get_range(nii, isLabel)
if size(nii.img, 8)>2 || any(nii.hdr.datatype == [128 511 2304]) % RGB / RGBA
    if max(nii.img(:))>2, rg = [0 255]; else, rg = [0 1]; end
    return;
elseif nii.hdr.cal_max~=0 && nii.hdr.cal_max>min(nii.img(:))
    rg = [nii.hdr.cal_min nii.hdr.cal_max];
    return;
end

img = nii.img(:,:,:,1);
img = img(:);
img(isnan(img) | isinf(img)) = [];
if ~isreal(img), img = abs(img); end
if ~isfloat(img)
    slope = nii.hdr.scl_slope; if slope==0, slope = 1; end
    img = single(img) * slope + nii.hdr.scl_inter;
end

mi = min(img); ma = max(img);
if nii.hdr.intent_code > 1000 || (nargin>1 && isLabel)
    rg = [mi ma]; return;
end

ind = abs(img)>50;
if sum(ind)<numel(img)/10, ind = abs(img)>std(img)/2; end
im = img(ind);
mu = mean(im);
sd = std(im);
rg = mu + [-2 2]*sd;
if rg(1)<=0, rg(1) = sd/5; end
if rg(1)<mi || isnan(rg(1)), rg(1) = mi; end
if rg(2)>ma || isnan(rg(2)), rg(2) = ma; end
if rg(1)==rg(2), rg(1) = mi; if rg(1)==rg(2), rg(1) = 0; end; end
% rg = round(rg, 2, 'significant'); % since 2014b
rg = str2num(sprintf('%.2g ', rg)); %#ok<*ST2NM>
if rg(1)==rg(2), rg(1) = mi; end
if abs(rg(1))>10, rg(1) = floor(rg(1)/2)*2; end % even number
if abs(rg(2))>10, rg(2) = ceil(rg(2)/2)*2; end % even number
end

