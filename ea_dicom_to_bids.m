function anat_files = ea_dicom_to_bids(subjID, fnames, dataset_folder)

% options that should appear in the table
table_options = {'skip','skip', 'skip';
    'anat','preop', 'T1w';
    'anat','preop', 'T2w';
    'anat','preop', 'FGATIR';
    'anat','preop', 'PDw';
    'anat','postop', 'CT';
    'anat','postop', 'ax_MR';
    'anat','postop', 'sag_MR';
    'anat','postop', 'cor_MR'};

% lookup table to pre-allocate options for the table.
lookup_table = loadjson(fullfile(ea_getearoot(), 'helpers', 'dicom_bids_lookuptable.json'));

nii_folder = fullfile(dataset_folder, 'sourcedata', subjID, 'tmp');    % where are the nifti files located?

imgs = cell(length(fnames),1);
h_wait = waitbar(0, 'Please wait while Niftii images are being loaded');
for image_idx = 1:length(fnames)
    imgs{image_idx} = struct();
    [imgs{image_idx}.p, imgs{image_idx}.frm, imgs{image_idx}.rg, imgs{image_idx}.dim] = read_nii(fullfile(nii_folder, [fnames{image_idx}, '.nii.gz']), [], 0);
    waitbar(image_idx / length(fnames), h_wait, sprintf('Please wait while Niftii images are being loaded (%i/%i)', image_idx, length(fnames)));
end
close(h_wait);

Session = categorical(repmat({'skip'},[length(fnames),1]),unique(table_options(:,2)));
Modality = categorical(repmat({'skip'},[length(fnames),1]), unique(table_options(:,3)));
Type = categorical(repmat({'anat'},[length(fnames),1]),unique(table_options(:,1)));
T = table(fnames, Session, Type, Modality);

try
T_preallocated = preallocate_table(T, lookup_table);
catch
disp('Preallocation of table failed, defaulting to skip!');
T_preallocated = T;
end

% create GUI
ui = dicom_to_bids;

% populate table
ui.niiFileTable.Data = T_preallocated;
ui.niiFileTable.ColumnEditable = [false true true, true];

% set subject ID and file path
ui.SubjectIDLabel.Text = subjID;
ui.FilepathLabel.Text = dataset_folder;

% update preview tree and expand it
ui.previewtree_subj.Text = subjID;
expand(ui.Tree, 'all');

preview_nii(ui, imgs{1,1}); % set initial image to the first one
update_preview_tree(ui, table_options, subjID) % call preview tree updater to get preallocated changes

ui.niiFileTable.CellSelectionCallback = @(src,event) preview_nii(ui,imgs{event.Indices(1), 1}); % callback for table selection -> display current selected image
ui.niiFileTable.CellEditCallback = @(src,event) update_preview_tree(ui, table_options, subjID); % callback for cell change -> update ui tree on the right

ui.UIFigure.WindowScrollWheelFcn = @(src, event) scroll_nii(ui, event);     % callback for scrolling images

% OK button behaviour
ui.OKButton.ButtonPushedFcn = @(btn,event) ok_button_function(ui, ui.niiFileTable, table_options, dataset_folder, nii_folder, subjID);

% cancel button behaviour
ui.CancelButton.ButtonPushedFcn =  @(btn,event) cancel_button_function(ui);

waitfor(ui.UIFigure);

try
    anat_files = getappdata(groot, 'anat_files');
catch
    anat_files = [];
    
end

end

function update_preview_tree(ui, table_options, subjID)

% update tree view as well
dat = cellfun(@char,table2cell(ui.niiFileTable.Data),'uni',0);  % get data
ui.previewtree_preop_anat.Children.delete;      % delete children
ui.previewtree_postop_anat.Children.delete;    % delete children

% populate tree
sessions = unique(table_options(:,2));
sessions(ismember(sessions,'skip')) = [];   % remove skip

for sesIdx = 1:length(sessions)
    
    ses = sessions{sesIdx};
    
    for fileIdx = 1:length(dat)
        if strcmp(dat{fileIdx, 2}, ses) && ~strcmp(dat{fileIdx, 4}, 'skip') && ~strcmp(dat{fileIdx, 3}, 'skip')
            
            fname = sprintf('%s_ses-%s_%s', subjID, ses, dat{fileIdx, 4});   % generate BIDS filename
            
            if strcmp(ses, 'preop')
                uitreenode(ui.previewtree_preop_anat, 'Text', fname);
            else
                uitreenode(ui.previewtree_postop_anat, 'Text', fname);
            end
        end
        
    end
end


end


function table_preallocated = preallocate_table(table, lookup_table)

table_preallocated = table;

% get sessions
image_types = fieldnames(lookup_table);

% filenames
for rowIdx = 1:height(table)
    
    fname = table.fnames{rowIdx};
    
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
                    end
                end
            end
        end
    end
end

end

function cancel_button_function(uiapp)

 s = uiconfirm(uiapp.UIFigure, 'Do you really want to cancel file selection?', 'Confirm close', ...
    'Options', {'Yes', 'No'}, 'Icon', 'question');

switch s
    case 'Yes'
        delete(uiapp);
end

    
end

function ok_button_function(uiapp, TT, table_options, dataset_folder, nii_folder, subjID)

dat = cellfun(@char,table2cell(TT.Data),'uni',0);
sessions = unique(table_options(:,2));
sessions(ismember(sessions,'skip')) = [];   % remove skip

% preallocate anat_files
anat_files = struct();
for sesIdx = 1:length(sessions)
    ses = sessions{sesIdx};
    anat_files.(ses) = [];
end

% before we do anything, first a sanity check if user has selected only one file per modality
sanity_check_passed = true;
for sesIdx = 1:length(sessions)
    
    ses = sessions{sesIdx};
    
    for fileIdx = 1:length(dat)
        if strcmp(dat{fileIdx, 2}, ses) && ~strcmp(dat{fileIdx, 4}, 'skip') && ~strcmp(dat{fileIdx, 3}, 'skip')
            
            fname = sprintf('%s_ses-%s_%s', subjID, ses, dat{fileIdx, 4});      % generate BIDS filename
            
            if ~isfield(anat_files.(ses), dat{fileIdx, 4})
                anat_files.(ses).(dat{fileIdx, 4}) = fname;                     % set output struct
            else % more than one file has been selected
                sanity_check_passed = false;
                s = uiconfirm(uiapp.UIFigure, 'Multiple files have been detected for one or more modalities. Please select only one file per modality and session', ...
                    'Too many files detected', 'Options', {'OK'}, ...
                    'Icon', 'warning');
                break
            end
            
        end
    end
    
end

% check if all images have been skipped
if all(any(ismember(dat(:,2:4),'skip'),2))
    s = uiconfirm(uiapp.UIFigure, 'All images have been skipped. Please select at least one preop and one postop image.', ...
        'No files selected', 'Options', {'OK'}, ...
        'Icon', 'warning');
    sanity_check_passed = false;
end

if sanity_check_passed == true
    for sesIdx = 1:length(sessions)
        
        ses = sessions{sesIdx};
        
        for fileIdx = 1:length(dat)
            if strcmp(dat{fileIdx, 2}, ses) && ~strcmp(dat{fileIdx, 4}, 'skip') && ~strcmp(dat{fileIdx, 3}, 'skip')
                
                fname = sprintf('%s_ses-%s_%s', subjID, ses, dat{fileIdx, 4});   % generate BIDS filename
                
                % move files
                if ~exist(fullfile(dataset_folder, 'rawdata', subjID, ['ses-', ses], 'anat'), 'dir')
                    mkdir(fullfile(dataset_folder, 'rawdata', subjID, ['ses-', ses], 'anat'));
                end
                copyfile(fullfile(nii_folder, [dat{fileIdx, 1}, '.nii.gz']), fullfile(dataset_folder, 'rawdata', subjID, ['ses-', ses], 'anat', [fname, '.nii.gz']));
                copyfile(fullfile(nii_folder, [dat{fileIdx, 1}, '.json']), fullfile(dataset_folder, 'rawdata', subjID, ['ses-', ses], 'anat', [fname, '.json']));
                
                anat_files.(ses).(dat{fileIdx, 4}) = fname; % set output struct
            end
            
        end
    end
end

if sanity_check_passed == true
    setappdata(groot, 'anat_files', anat_files);
    delete(uiapp);      % close window
end

end


function preview_nii(ui, img)

% update info area
try
    time_and_date_ugly = img.p.nii.ext.edata_decoded.AcquisitionDateTime;
    time_and_date_pretty = sprintf('%s.%s.%s %s:%s', num2str(time_and_date_ugly(7:8)), ...
        num2str(time_and_date_ugly(5:6)), num2str(time_and_date_ugly(1:4)), ...
        num2str(time_and_date_ugly(9:10)), num2str(time_and_date_ugly(11:12)));
catch
    time_and_date_pretty = 'N/A';
end
info_str = sprintf('Size: [%s x %s x %s]\nPixel dimensions: [%.2f x %.2f x %.2f]\nAcquistion date: %s\nIntensity range: [%.0f, %.0f]', ...
    num2str(img.dim(1)), num2str(img.dim(2)), num2str(img.dim(3)), ...
    img.p.pixdim(1), img.p.pixdim(2), img.p.pixdim(3), ...
    time_and_date_pretty, ...
    min(img.p.nii.img(:)), max(img.p.nii.img(:)));

ui.infoArea.Value = {info_str};

% plot images
setappdata(ui.UIFigure, 'img', img);

% axial
cut_slice = round(img.dim(3)/2);
imagesc(ui.axes_axi, img.p.nii.img(:, :, cut_slice));
ui.axes_axi.Colormap = gray(128);
setappdata(ui.UIFigure, 'cut_slice_axi', cut_slice); % save current cut slice for scrolling
ui.axes_axi.DataAspectRatioMode = 'manual';
ui.axes_axi.DataAspectRatio = [img.p.pixdim(1), img.p.pixdim(2), 1];
set(ui.axes_axi, 'view', [90, -90]);

% coronal
cut_slice = round(img.dim(2)/2);
imagesc(ui.axes_cor, squeeze(img.p.nii.img(:, cut_slice, :)));
ui.axes_cor.Colormap = gray(128);
setappdata(ui.UIFigure, 'cut_slice_cor', cut_slice); % save current cut slice for scrolling
ui.axes_cor.DataAspectRatioMode = 'manual';
ui.axes_cor.DataAspectRatio = [img.p.pixdim(1), img.p.pixdim(3), 1];
set(ui.axes_cor, 'view', [90, -90]);

% sagittal
cut_slice = round(img.dim(1)/2);
imagesc(ui.axes_sag, squeeze(img.p.nii.img(cut_slice, :, :)));
ui.axes_sag.Colormap = gray(128);
setappdata(ui.UIFigure, 'cut_slice_sag', cut_slice); % save current cut slice for scrolling
ui.axes_sag.DataAspectRatioMode = 'manual';
ui.axes_sag.DataAspectRatio = [img.p.pixdim(1), img.p.pixdim(3), 1];
set(ui.axes_sag, 'view', [90, -90]);

end

function scroll_nii(ui, event)

hAxes = checkMousePointer(ui.UIFigure, ui.CenterPanel);
img = getappdata(ui.UIFigure, 'img');
dim = img.dim;

if ~isempty(hAxes)
    switch hAxes.Tag
        case 'axi'
            sliceNr = getappdata(ui.UIFigure, 'cut_slice_axi');
            if event.VerticalScrollCount == -1 % up scroll
                if ~(sliceNr >= img.dim(3) - 2)
                    sliceNr = sliceNr + 2;
                end
            else
                if ~(sliceNr <= 2)
                    sliceNr = sliceNr - 2;
                end
                
            end
            imagesc(ui.axes_axi, img.p.nii.img(:, :, sliceNr));
            setappdata(ui.UIFigure, 'cut_slice_axi', sliceNr);
        case 'cor'
            sliceNr = getappdata(ui.UIFigure, 'cut_slice_cor');
            if event.VerticalScrollCount == -1 % up scroll
                if ~(sliceNr >= img.dim(2) - 2)
                    sliceNr = sliceNr + 2;
                end
            else
                if ~(sliceNr <= 2)
                    sliceNr = sliceNr - 2;
                end
                
            end
            imagesc(ui.axes_cor, squeeze(img.p.nii.img(:, sliceNr, :)));
            setappdata(ui.UIFigure, 'cut_slice_cor', sliceNr);
        case 'sag'
            sliceNr = getappdata(ui.UIFigure, 'cut_slice_sag');
            if event.VerticalScrollCount == -1 % up scroll
                if ~(sliceNr >= img.dim(1) - 2)
                    sliceNr = sliceNr + 2;
                end
            else
                if ~(sliceNr <= 2)
                    sliceNr = sliceNr - 2;
                end
                
            end
            imagesc(ui.axes_sag, squeeze(img.p.nii.img(sliceNr, :, :)));
            setappdata(ui.UIFigure, 'cut_slice_sag', sliceNr);
        otherwise
            
    end
end

end

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

