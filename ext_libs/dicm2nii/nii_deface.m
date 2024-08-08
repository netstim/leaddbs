function nii_deface(in)
% Syntax: NII_DEFACE(in)
%  input can be a single anat 3D NIfTI, or BIDS folder containing sub-* folders.
%  If no input is provided, a file/folder dialog will pop up.
% 
% NII_DEFACE removes face and neck structure from T1w/T2w NIfTI.
% 
% How does this code deface? 
%  1. A xform matrix is estimated for normalizing an anat image to MNI template
%  2. MNI deface mask, which includes brain and scalp but without neck and face,
%     is transformed to the subject reference using the above xform matrix
%  3. Voxels outside the transformed mask is set to zeros
%  4. A GUI gives options to overwrite the file etc.
%  
% Potential benefit over other defacing tools:
%  1. Except zeroing some voxels, nothing else of the NIfTI is altered.
%  2. Most neck tissue is removed, which reduces the chance of error for the 
%     anatomical analysis for some tools.
%  3. File size is significantly reduced for .gz version. 
%
% See also NII_TOOL NII_COREG NII_VIEWER NII_XFORM

% 230415 Wrote it by Xiangrui.Li@gmail.com

global DEFACE_PARAMS %#ok<*GVMIS>
DEFACE_PARAMS = struct('noui', false, 'tag', '', 'dir', '');

if nargin<1 || isempty(in)
    jFileChooser = dicm2nii('', 'jFileChooser', 'func_handle');
    in = jFileChooser([], 'Select a dataset folder or a NIfTI file', false);
    if isnumeric(in), return; end
end

f = fileparts(which('nii_viewer'));
niiT1 = nii_tool('load', [f '/templates/MNI_2mm_T1.nii']);
niiT2 = nii_tool('load', [f '/templates/MNI_2mm_T2.nii']);
msk = load([f '/example_data.mat'], 'MNI_2mm_deface_mask');
msk = msk.MNI_2mm_deface_mask; % enlarged brain+scalp mask

if isfile(in)
    deface1(in, niiT1, niiT2, msk);
    return;
end

f = [in filesep];
subjs = dir([f 'sub-*']); % BIDS folder structure
for i = 1:numel(subjs)
    subj = subjs(i).name;
    nams = dir([f subj '/anat/' subj '*.nii.gz']);    
    for j = 1:numel(nams)
        fnam = [nams(j).folder filesep nams(j).name];
        deface1(fnam, niiT1, niiT2, msk);
    end
end

%% Do one brain
function deface1(fnam, niiT, niiT2, msk)
global DEFACE_PARAMS
niiM = nii_tool('load', fnam);
if endsWith(fnam, 'T2w.nii.gz'), niiT = niiT2; end
[M, mss] = nii_coreg(niiT, niiM);
if mss>0.5, warning('nii_deface:AlignBad', 'Alignment unreliable: %s', fnam); end
M = M \ [msk.hdr.sform_mat; 0 0 0 1];
msk = nii_tool('update', msk, M);
msk.hdr.sform_code = niiM.hdr.sform_code;
msk = nii_xform(msk, niiM.hdr, [], 'nearest', 0);
nii = niiM;
slp = nii.hdr.scl_slope; if slp==0, slp = 1; end
nii.img(~msk.img) = -nii.hdr.scl_inter/slp; % set outside to 0

% how to save the result?
if DEFACE_PARAMS.tag == "another"
    subj = regexp(fnam, '(?<=.*[/\\])sub-\w+(?=[/\\]anat)', 'match', 'once');
    s  = ['\' filesep];
    DEFACE_PARAMS.dir = regexprep(DEFACE_PARAMS.dir, ...
        '[/\\]sub-\w+[/\\]anat', [s subj s 'anat'], 'once');
end

if ~DEFACE_PARAMS.noui
    fh = nii_viewer(niiM, nii); % overlay for check
    waitfor(dialog_show(fh.Position));
    try close(fh); catch, end
    if DEFACE_PARAMS.tag=="another"
        pth = uigetdir(DEFACE_PARAMS.dir, 'Select folder to save the defaced NIfTI');
        if ~isnumeric(pth), DEFACE_PARAMS.dir = pth; end
    end
end

if DEFACE_PARAMS.tag == "overwrite"
    nii_tool('save', nii);
elseif DEFACE_PARAMS.tag ==  "another"
    [~, nam, ext] = fileparts(niiM.hdr.file_name);
    nam = [DEFACE_PARAMS.dir filesep nam ext];
    if ~isfolder(DEFACE_PARAMS.dir), mkdir(DEFACE_PARAMS.dir); end
    nii_tool('save', nii, nam);
elseif DEFACE_PARAMS.tag ==  "bad"
    fprintf(2, '%s\n', ['Skipped: ' fnam]);
elseif DEFACE_PARAMS.noui
    fh = nii_viewer(niiM, nii);
    uiwait(fh, 20); % stay until user close nii_viewer or 20 seconds
    try close(fh); catch, end
end

%% Action dialog
function d = dialog_show(pos)
d = dialog('Name', 'Save the defaced?', 'CloseRequestFcn', @dialog_cb, ...
    'Position', [pos(1)+pos(3) pos(2) 280 180], 'WindowStyle', 'normal');
uicontrol(d, 'Style', 'checkbox', 'Position', [40 140 240 24], 'Value', 0, ...
    'String', 'Do the same for remaining subjects');
uiBtn = @(y,tag,str)uicontrol(d, 'Position', [40 y 180 24], 'Tag', tag, ...
    'String', str, 'Callback', @dialog_cb);
uiBtn(20, 'overwrite', 'Overwrite the source file');
uiBtn(50, 'another', 'Save to another BIDS folder');
uiBtn(80, 'nosave', 'No save, just check one by one');
uiBtn(110, 'bad', 'Bad result & disp file name');

%% dialog callback
function dialog_cb(h, evt)
global DEFACE_PARAMS
fh = ancestor(h, 'figure');
ck = findobj(fh, 'Style', 'checkbox');
DEFACE_PARAMS.noui = ck.Value;
if evt.EventName == "Close"
    DEFACE_PARAMS.tag = '';
else
    DEFACE_PARAMS.tag = h.Tag;
    if h.Tag == "bad", DEFACE_PARAMS.noui = false; end
end
delete(fh);

%%