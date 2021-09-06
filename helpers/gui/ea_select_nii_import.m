function ea_select_nii_import(fnames, nii_folder, rawdata_folder)

valueset = {'skip','skip', 'skip';
    'anat','preop', 'T1w';
    'anat','preop', 'T2w';
    'anat','preop', 'FGATIR';
    'anat','preop', 'PDw';
    'anat','postop', 'CT';
    'anat','postop', 'ax_MR';
    'anat','postop', 'sag_MR';
    'anat','postop', 'cor_MR'};

verLessThanOctave = false; % debug

imgs = cell(length(fnames),1);
h_wait = waitbar(0, 'Please wait while Niftii images are being loaded ');
for image_idx = 1:length(fnames)
    imgs{image_idx} = struct();
    [imgs{image_idx}.p, imgs{image_idx}.frm, imgs{image_idx}.rg, imgs{image_idx}.dim] = read_nii(fullfile(nii_folder, [fnames{image_idx}, '.gz']));
    waitbar(image_idx / length(fnames));
end
close(h_wait);

% debug
Subject = {'01'};
Session = {'01'};
AcquisitionDate = {'N/A'};
Comment = {'N/A'};
S = table(Subject,Session,AcquisitionDate,Comment);

Session = categorical(repmat({'skip'},[length(fnames),1]),unique(valueset(:,2)));
Modality = categorical(repmat({'skip'},[length(fnames),1]), unique(valueset(:,3)));
Type = categorical(repmat({'skip'},[length(fnames),1]),unique(valueset(:,1)));
T = table(fnames, Session, Type, Modality);

ModalityTablePref = getpref('dicm2nii_gui_para', 'ModalityTable', T);

% GUI
figsize = [0.8, 0.6];

setappdata(0,'Canceldicm2nii',false)
scrSz = get(0, 'ScreenSize');
mpos = get(0, 'MonitorPositions');
clr = [1 1 1]*206/256;
figargs = {'Units', 'normalized', ...
    'Position',[0.1 0.2 figsize(1) figsize(2)],...
    'Color', clr,...
    'CloseRequestFcn',@my_closereq};
if verLessThanOctave
    hf = figure(figargs{1});
    set(hf,figargs{2:end});
    % add help
    set(hf,'ToolBar','none')
    set(hf,'MenuBar','none')
else
    hf = uifigure(figargs{:});
end
uimenu(hf,'Label','help','Callback',@(src,evnt) showHelp(valueset))
set(hf,'Name', 'dicm2nii - BIDS Converter', 'NumberTitle', 'off')

% tables
if verLessThanOctave
    SCN = S.Properties.VariableNames;
    S   = table2cell(S);
    TCN = T.Properties.VariableNames;
    T   = cellfun(@char,table2cell(T),'uni',0);
end
TS = uitable(hf,'Data',S);
set(TS, 'Units', 'normalized');
TT = uitable(hf,'Data',T);
set(TT, 'Units', 'normalized');

TSpos = [0 0.8 0.6 0.2];
TTpos = [0 0.1 0.6 0.6];

TS.Position = TSpos;
TT.Position = TTpos;

TS.ColumnEditable = [true true true true];

TT.ColumnEditable = [false true true, true];
setappdata(0,'ModalityTable',TT.Data)
setappdata(0,'SubjectTable',TS.Data)

% button
Bpos = [hf.Position(3)-120 20 100 30];
BCB  = @(btn,event) BtnModalityTable(hf,TT, TS);

B = uibutton(hf,'Position',Bpos);
B.Text = 'OK';
B.ButtonPushedFcn = BCB;


% preview panel
axesArgs_axi = {hf,'Units', 'normalized', 'Position',[0.55 0.5 0.4 0.4],...
    'Colormap',gray(64)};
ax_axi = preview_nii([],imgs{1, 1}, axesArgs_axi, 'axi', hf);
ax_axi.YTickLabel = [];
ax_axi.XTickLabel = [];

% coronal slice
%axesArgs_cor = {hf,'Position',[hf.Position(3) 100 100 hf.Position(4)-90],...
%    'Colormap',gray(64)};
%ax_cor = preview_nii([],imgs{1, 1}, axesArgs_cor, 'cor');
%ax_cor.YTickLabel = [];
%ax_cor.XTickLabel = [];

% sagittal slice
%axesArgs_sag = {hf,'Position',[hf.Position(3)-120 120 100 hf.Position(4)-90],...
%    'Colormap',gray(64)};
%ax_sag = preview_nii([],imgs{1, 1}, axesArgs_sag, 'sag');
%ax_sag.YTickLabel = [];
%ax_sag.XTickLabel = [];

TT.CellSelectionCallback = @(src,event) preview_nii(ax_axi,imgs{event.Indices(1), 1}, axesArgs_axi, 'axi', hf);
hf.WindowScrollWheelFcn = @(src, event) scroll_nii(ax_axi, event);  % this callback enables scrolling

waitfor(hf);
if getappdata(0,'Canceldicm2nii')
    return;
end

% get results
ModalityTable = getappdata(0,'ModalityTable');
SubjectTable = getappdata(0,'SubjectTable');

% setpref
if istable(ModalityTable)
    ModalityTable = cellfun(@char,table2cell(ModalityTable),'uni',0);
end
ModalityTableSavePref = ModalityTable(~any(ismember(ModalityTable(:,2:3),'skip'),2),:);

for imod = 1:size(ModalityTableSavePref,1)
    match = cellfun(@(Mod) strcmp(Mod,ModalityTableSavePref{imod,1}),table2cell(ModalityTablePref(:,1)));
    if any(match) % replace old pref
        ModalityTablePref.Type(match) = ModalityTableSavePref{imod,2};
        ModalityTablePref.Modality(match) = ModalityTableSavePref{imod,3};
    else % append new pref
        ModalityTablePref = [ModalityTablePref;ModalityTableSavePref(imod,:)];
    end
    
    setpref('dicm2nii_gui_para', 'ModalityTable', ModalityTablePref);
    
end
end

function scroll_nii(ax, event)
    
    img = getappdata(ax, 'img');
    sliceNr = getappdata(ax, 'cut_slice');
    
   if event.VerticalScrollCount == -1 % up scroll
       sliceNr = sliceNr + 1;
   else
       sliceNr = sliceNr - 1;
   end
    imagesc(ax, img.p.nii.img(:, :, sliceNr));
    setappdata(ax, 'cut_slice', sliceNr);
    
end

function ax = preview_nii(ax, img, axesArgs, cut_direction, dicom_fig)

if isempty(ax)
    ax=uiaxes(axesArgs{:});
end

if isempty(getappdata(dicom_fig, 'infoBox'))
    infoBox = annotation(dicom_fig, 'textbox', [0.55 0.2 0.4 0.4], 'String', ' ');
    setappdata(dicom_fig, 'infoBox', infoBox);
else
    infoBox = getappdata(dicom_fig, 'infoBox');
end

% create info panel
try
time_and_date_ugly = img.p.nii.ext.edata_decoded.AcquisitionDateTime;
time_and_date_pretty = sprintf('%s.%s.%s %s:%s', num2str(time_and_date_ugly(7:8)), ...
    num2str(time_and_date_ugly(5:6)), num2str(time_and_date_ugly(1:4)), ...
    num2str(time_and_date_ugly(9:10)), num2str(time_and_date_ugly(11:12)));
catch
    time_and_date_pretty = 'N/A';
end
info_str = sprintf('Size: [%s x %s %s]\nPixel dimensions: [%.2f x %.2f x %.2f]\nAcquistion Date: %s\n', ... 
    num2str(img.dim(1)), num2str(img.dim(2)), num2str(img.dim(3)), ...
    img.p.pixdim(1), img.p.pixdim(2), img.p.pixdim(3), ...
    time_and_date_pretty);

infoBox.String = info_str;

setappdata(ax, 'img', img);     % save .nii image in current appdata ax for scrolling

axis(ax,'off');
set(ax,'BackgroundColor',[0 0 0])

if strcmp('axi', cut_direction)
    cut_slice = round(img.dim(3)/2);
    imagesc(ax, img.p.nii.img(:, :, cut_slice));
elseif strcmp('sag', cut_direction)
    cut_slice = round(img.dim(1)/2);
    imagesc(ax, permute(img.p.nii.img(cut_slice, :, :), [2 3 1]));
else
    cut_slice = round(img.dim(2)/2);
    imagesc(ax, permute(img.p.nii.img(:, cut_slice, :), [1 3 2]));
end
ax.DataAspectRatio = [img.p.pixdim(1), img.p.pixdim(2), 1];

setappdata(ax, 'cut_slice', cut_slice); % save current cut slice for scrolling

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

%% reorient nii to diagnal major
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

%% Estimate lower and upper bound of img display
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

%%
function BtnModalityTable(h,TT,TS)
dat = cellfun(@char,table2cell(TT.Data),'uni',0);

if all(any(ismember(dat(:,2:3),'skip'),2))
    warndlg('All images are skipped... Please select the type and modality for all scans','No scan selected');
    return;
end
setappdata(0,'ModalityTable',TT.Data)
setappdata(0,'SubjectTable',TS.Data)
delete(h)
end

function my_closereq(src,~)
% Close request function
% to display a question dialog box

selection = uiconfirm(src,'Cancel Dicom conversion?',...
    'Close dicm2nii');

switch selection
    case 'OK'
        delete(src)
        setappdata(0,'Canceldicm2nii',true)
    case 'Cancel'
        return
end
end



