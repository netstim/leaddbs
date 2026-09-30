function app = ea_prepare_legui(app) %#ok<INUSD>
%preferBrainshiftCT = false;
enforceSameGrid    = false;
useMNI             = false;

% ------- Determine subject root (Lead-DBS subject folder) -------
%     app.SelDir = LeG_lastDir();
%     subjRoot = app.SelDir;
%     [~, subjId] = fileparts(subjRoot);

subjRoot = app.SelDir;

% Sometimes the folder path is stored inside a one-item cell
if iscell(subjRoot) && numel(subjRoot) == 1
    subjRoot = subjRoot{1};
end

% Make sure it is normal text
subjRoot = char(subjRoot);
[~, subjId] = fileparts(subjRoot);

% ------- Determine whether approved brain shift should be used -------
brainshiftApproved = isBrainshiftApproved_local(subjRoot, subjId);

if brainshiftApproved
    fprintf('Brain shift is APPROVED for %s.\n', subjId);
else
    fprintf('Brain shift is not approved or was not run for %s.\n', subjId);
end

GM_filename = strcat(subjId, '_ses-preop_space-anchorNative_desc-preproc_acq-iso_label-GM_mod-iso_T1w_mask.nii');
WM_filename = strcat(subjId, '_ses-preop_space-anchorNative_desc-preproc_acq-iso_label-WM_mod-iso_T1w_mask.nii');
CSF_filename = strcat(subjId, '_ses-preop_space-anchorNative_desc-preproc_acq-iso_label-CSF_mod-iso_T1w_mask.nii');
Surfaces_filename = strcat(subjId, '_ses-preop_space-anchorNative_desc-preproc_acq-iso_surface.mat');
Electrodes_filename = strcat(subjId, '_electrodes.mat');
Channels_filename = strcat(subjId, '_channels.mat');
% ------- Find MR/CT inside that Lead-DBS tree -------
try
    [mrFile, ctFile] = findLeadDBSImagePair_local(subjRoot, useMNI);
catch ME
    msgbox(sprintf('Could not find MR/CT in:\n%s\n\n%s', subjRoot, ME.message));
    app.LoadImgsBtnH.Enable = "on";
    return;
end

app.DICOMDir = '';  % unused downstream
[~, app.SelFolder] = fileparts(app.SelDir);
% ------- UI label -------
try
    [~,mrBase] = fileparts(mrFile);
    app.PatientIDStr    = sprintf('%s (Lead-DBS)', mrBase);
    app.MainFigure.Name = sprintf('LeGUI v%0.1f - %s', app.Version, app.PatientIDStr);
catch
end

% ------- Load MR/CT and prepare GM/WM segmentations -------
app.WaitH = uiprogressdlg(app.MainFigure,'Indeterminate','on','Message','Loading MR/CT...');
app.MRInfo = spm_vol(mrFile);
app.MRImg = spm_read_vols(app.MRInfo);
app.CTInfo = spm_vol(ctFile);
app.CTImg = spm_read_vols(app.CTInfo);

%     [MRInfo, MRImg, CTInfo, CTImg] = loadReslicedMRCT(mrFile, ctFile);
%     app.MRInfo = MRInfo;
%     app.MRImg  = MRImg;
%     app.CTInfo = CTInfo;
%     app.CTImg  = CTImg;

[app.MRImg, MRRng] = normalizeImage(app.MRImg);
[app.CTImg, app.CTRng] = normalizeImage(app.CTImg);

if exist(fullfile(subjRoot, 'coregistration', 'anat', GM_filename), 'file') && ...
        exist(fullfile(subjRoot, 'coregistration', 'anat', WM_filename), 'file') && ...
        exist(fullfile(subjRoot, 'coregistration', 'anat', CSF_filename), 'file')

    % load existing GM/WM
    GM = spm_vol(char(fullfile(subjRoot, 'coregistration', 'anat', GM_filename)));
    WM = spm_vol(char(fullfile(subjRoot, 'coregistration', 'anat', WM_filename)));
    CSF = spm_vol(char(fullfile(subjRoot, 'coregistration', 'anat', CSF_filename)));

else
    % run segmentation
    [GM, WM, CSF] = get_segmentations(mrFile, subjRoot, GM_filename, WM_filename, CSF_filename);
end
app.GrayInfo = GM;
app.WhiteInfo = WM;
app.CSFInfo = CSF;
app.GrayImg = spm_read_vols(app.GrayInfo); 
app.WhiteImg = spm_read_vols(app.WhiteInfo);
app.CSFImg = spm_read_vols(app.CSFInfo);
SR = [];
app.BrainSurfRaw = []; app.ProjSurfRaw = [];
[app.BrainSurfRaw,app.ProjSurfRaw] = LeG_genSurfaces({app.GrayImg,app.WhiteImg},app.GrayInfo);
SR.BrainSurfRaw = app.BrainSurfRaw;
SR.ProjSurfRaw = app.ProjSurfRaw;
ea_mkdir(fullfile(subjRoot, 'surfaces'));
app.SurfacesFile = fullfile(subjRoot, 'surfaces', Surfaces_filename);
save(app.SurfacesFile,'-struct','SR')

% ------- Atlas / scaling / initial display -------
app.XYZScale   = sqrt(sum(app.MRInfo.mat(1:3,1:3).^2));
app.CurSlice   = round(size(app.MRImg, app.CurView2DIdx)/2);
app.DispImg    = app.MRImg;
app.DispImgSub = app.CTImg;
set(app.TotalSliceTextH,'text', ['/ ' num2str(size(app.DispImg, app.CurView2DIdx))]);

%---------- Store normalization matrices within app ------------
transDir = fullfile(subjRoot, 'normalization', 'transformations');
% look for files containing the convention
mni2nat = dir(fullfile(transDir, '*MNI152NLin2009bAsym_to-anchorNative*'));
nat2mni = dir(fullfile(transDir, '*anchorNative_to-MNI152NLin2009bAsym*'));

if ~isempty(mni2nat)
    app.yMR = fullfile(transDir, mni2nat(1).name);
    app.iyMR = fullfile(transDir, nat2mni(1).name);
else
    warning('No matching transformation file found.');
    app.yMR = '';
end

%------------ Saving filepath for electrodes and channelmap files -----------------
app.Reco = fullfile(subjRoot, 'reconstruction');
app.subjId = strcat(subjId, '_desc-reconstruction.mat');
app.recofile = fullfile(app.Reco, Electrodes_filename); %'Electrodes.mat'); %fullfile(app.Reco, app.subjId);
app.channelfile = fullfile(app.Reco, Channels_filename); %'Channel.mat');

% ------- Electrodes / ChannelMap -------
app.WaitH.Message = 'Loading electrodes...';
app.ElecImgProjRaw = zeros(size(app.MRImg), 'uint8');

% Keep both historical and current property names synchronized
app.ElectFile = app.recofile;
app.ChanMapFile = app.channelfile;


if isfile(app.ElectFile)

    ES = load(app.ElectFile);

    copyIf_local(app, ES, { ...
        'ElecXYZRaw', ...
        'ElecFullIdxRaw', ...
        'ElecCOMIdxRaw', ...
        'ElecXYZProjRaw', ...
        'ElecFullIdxProjRaw', ...
        'ElecCOMIdxProjRaw', ...
        'ElecMapRaw', ...
        'DepthElecRaw', ...
        'MicroElecRaw', ...
        'RefElecRaw', ...
        'GndElecRaw', ...
        'ShaftNames', ...
        'ShaftModels', ...
        'ShaftMembership'});

    if isfield(ES,'PatientIDStr') && ...
            ~strcmp(app.PatientIDStr, ES.PatientIDStr)

        app.PatientIDStr = ES.PatientIDStr;

        app.MainFigure.Name = sprintf( ...
            'LeGUI v%0.1f - %s', ...
            app.Version, ...
            app.PatientIDStr);
    end

else
    warning('Electrodes file was not found: %s', app.ElectFile);
end

% restore contacts from Lead-DBS reconstruction.mat

reconstructionFile = fullfile(app.Reco, app.subjId);

if isempty(app.ElecXYZProjRaw) && isfile(reconstructionFile)

    R = load(reconstructionFile, 'reco');

    if isfield(R, 'reco') && ...
            isfield(R.reco, 'native') && ...
            isfield(R.reco.native, 'coords_mm') && ...
            ~isempty(R.reco.native.coords_mm)

        coordsCells = R.reco.native.coords_mm;
        coordsCells = coordsCells(~cellfun(@isempty, coordsCells));

        if ~isempty(coordsCells)

            % Combine all electrode shafts into one N-by-3 contact list
            allCoords = vertcat(coordsCells{:});

            % Restore LeGUI world-coordinate arrays
            app.ElecXYZRaw     = allCoords;
            app.ElecXYZProjRaw = allCoords;

            % Convert world coordinates to MR voxel indices
            voxelCoords = app.MRInfo.mat \ ...
                [allCoords'; ones(1, size(allCoords,1))];

            voxelCoords = round(voxelCoords(1:3,:)');

            app.ElecCOMIdxRaw     = voxelCoords;
            app.ElecCOMIdxProjRaw = voxelCoords;

            % Rebuild the voxel spheres used by LeGUI
            app.ElecFullIdxRaw     = cell(size(allCoords,1),1);
            app.ElecFullIdxProjRaw = cell(size(allCoords,1),1);

            radiusMM = str2double(app.ElecRadEditH.Value);

            if ~isfinite(radiusMM) || radiusMM <= 0
                radiusMM = 1.3;
            end

            searchRadius = ceil(radiusMM ./ app.XYZScale);
            imageSize = size(app.MRImg);

            for contactIndex = 1:size(voxelCoords,1)

                center = voxelCoords(contactIndex,:);

                xRange = max(1,center(1)-searchRadius(1)): ...
                    min(imageSize(1),center(1)+searchRadius(1));

                yRange = max(1,center(2)-searchRadius(2)): ...
                    min(imageSize(2),center(2)+searchRadius(2));

                zRange = max(1,center(3)-searchRadius(3)): ...
                    min(imageSize(3),center(3)+searchRadius(3));

                [X,Y,Z] = ndgrid(xRange,yRange,zRange);

                distanceMM = sqrt( ...
                    ((X-center(1))*app.XYZScale(1)).^2 + ...
                    ((Y-center(2))*app.XYZScale(2)).^2 + ...
                    ((Z-center(3))*app.XYZScale(3)).^2);

                insideSphere = distanceMM <= radiusMM;

                sphereIndices = [ ...
                    X(insideSphere), ...
                    Y(insideSphere), ...
                    Z(insideSphere)];

                app.ElecFullIdxRaw{contactIndex} = sphereIndices;
                app.ElecFullIdxProjRaw{contactIndex} = sphereIndices;
            end

            % Ensure the electrode assignment table has the correct rows
            if isempty(app.ElecMapRaw) || ...
                    size(app.ElecMapRaw,1) ~= size(allCoords,1)

                app.ElecMapRaw = [ ...
                    repmat({'NaN'},size(allCoords,1),1), ...
                    num2cell(nan(size(allCoords,1),2))];
            end

            fprintf( ...
                'Restored %d contacts from %s\n', ...
                size(allCoords,1), reconstructionFile);
        end
    end
end
% Load the saved Assign GUI tables from the channels file

if isfile(app.ChanMapFile)

    CM = load(app.ChanMapFile);

    if isfield(CM, 'LabelMap')
        app.LabelMap = CM.LabelMap;
    else
        warning('LabelMap is missing from: %s', app.ChanMapFile);
    end

    if isfield(CM, 'ChannelMap1')
        app.ChannelMap1 = CM.ChannelMap1;
    else
        warning('ChannelMap1 is missing from: %s', app.ChanMapFile);
    end

    if isfield(CM, 'ChannelMap2')
        app.ChannelMap2 = CM.ChannelMap2;
    else
        warning('ChannelMap2 is missing from: %s', app.ChanMapFile);
    end

else
    warning('Channels file was not found: %s', app.ChanMapFile);
end

fprintf(['Loaded %d contacts, %d assignment columns, ' ...
    '%d shaft names, and %d shaft models.\n'], ...
    size(app.ElecMapRaw,1), ...
    size(app.LabelMap,2), ...
    numel(app.ShaftNames), ...
    numel(app.ShaftModels));
end
% -------------------- helpers (local scope) --------------------

function [mrFile, ctFile] = findLeadDBSImagePair_local(rootDir, useMNI)
if useMNI
    mrPats = [ ...
        "normalization/anat/*_space-MNI152NLin2009bAsym_desc-preproc_acq-*_T1w.nii", ...
        "normalization/anat/*_space-MNI152NLin2009bAsym_desc-preproc_T1w.nii" ...
        ];
    ctPats = [ ...
        "normalization/anat/*_ses-postop_space-MNI152NLin2009bAsym_desc-preproc_CT.nii" ...
        ];
else
    mrPats = [ ...
        "coregistration/anat*/**/*_space-anchorNative_*_T1w.nii", ...
        "normalization/anat/*_acq-*_T1w.nii" ...
        ];
   
    ctPats = [ ...
        "coregistration/anat*/**/*_ses-postop_space-anchorNative_*desc-preproc_CT.nii", ...
        "preprocessing/*_ses-postop_space-anchorNative_desc-preproc_CT.nii" ...
        ];
end


mrFile = firstMatch_local(rootDir, mrPats);
ctFile = firstMatch_local(rootDir, ctPats);

assert(mrFile~="", 'Lead-DBS MR not found under: %s', rootDir);
assert(ctFile~="", 'Lead-DBS CT not found under: %s', rootDir);

mrFile = char(mrFile); ctFile = char(ctFile);
end

function p = firstMatch_local(rootDir, patterns)
p = "";
for k = 1:numel(patterns)
    hits = dir(fullfile(rootDir, patterns(k)));
    hits = hits(~[hits.isdir]);
    if ~isempty(hits)
        % pick the shortest path (usually canonical)
        [~, idx] = min(arrayfun(@(h) strlength(fullfile(h.folder,h.name)), hits));
        p = string(fullfile(hits(idx).folder, hits(idx).name));
        return;
    end
end
end

function copyIf_local(app, S, fields)
for i = 1:numel(fields)
    if isfield(S, fields{i}), app.(fields{i}) = S.(fields{i}); end
end
end

function [GM, WM, CSF] = get_segmentations(mrPath, subjRoot, GM_filename, WM_filename, CSF_filename)
%SEGMENT_GM_WM_INMEMORY Run SPM12 segmentation and return GM/WM masks in memory.
%
% Inputs
%   mrPath     : full path to MRI NIfTI to segment (already coregistered)
%   probThresh : probability threshold for binarizing c1/c2 (default 0.5)
%
% Outputs
%   GMmask : binary gray matter mask (logical 3D array)
%   WMmask : binary white matter mask (logical 3D array)
%   c1Vol  : GM probability map (double array)
%   c2Vol  : WM probability map (double array)
% --- SPM defaults ---
spm('defaults','fmri');
spm_jobman('initcfg');

% --- Build minimal segmentation struct ---
seg.channel.vols     = {mrPath};
seg.channel.biasreg  = 0.001;
seg.channel.biasfwhm = 60;
seg.channel.write    = [0 0]; % don't write bias-corrected MRI

tpm = fullfile(spm('Dir'),'tpm','TPM.nii');
ngausVals = [1 1 2 3 4 2];

for i = 1:6
    seg.tissue(i).tpm    = {sprintf('%s,%d', tpm, i)};
    seg.tissue(i).ngaus  = ngausVals(i);
    seg.tissue(i).native = [i<=3, 0]; % GM/WM/CSF probability maps
    seg.tissue(i).warped = [0 0];
end

seg.warp.mrf     = 1;
seg.warp.cleanup = 1;
seg.warp.reg     = [0 0.001 0.5 0.05 0.2];
seg.warp.affreg  = 'mni';
seg.warp.fwhm    = 0;
seg.warp.samp    = 3;
seg.warp.write   = [0 0];

% --- Run segmentation ---
fout = spm_preproc_run(seg);

% --- Load GM/WM volumes directly ---
c1Path = fout.tiss(1).c{1}; % GM prob map
c2Path = fout.tiss(2).c{1}; % WM prob map
c3Path = fout.tiss(3).c{1};
%     GM = spm_vol(c1Path); %c1Vol = spm_read_vols(V1);
%     WM = spm_vol(c2Path); %c2Vol = spm_read_vols(V2);
%     CSF = spm_vol(c3Path);
%     copyfile(c1Path, fullfile(subjRoot, 'coregistration', 'anat', GM_filename));
%     copyfile(c2Path, fullfile(subjRoot, 'coregistration', 'anat', WM_filename));
%     copyfile(c3Path, fullfile(subjRoot, 'coregistration', 'anat', CSF_filename));
%     ea_delete(c1Path);
%     ea_delete(c2Path);
%     ea_delete(c3Path);
% Permanent destination paths
gmDest = fullfile(subjRoot, 'coregistration', 'anat', GM_filename);
wmDest = fullfile(subjRoot, 'coregistration', 'anat', WM_filename);
csfDest = fullfile(subjRoot, 'coregistration', 'anat', CSF_filename);

% Copy temporary SPM outputs to the Lead-DBS subject folder
copyfile(c1Path, gmDest);
copyfile(c2Path, wmDest);
copyfile(c3Path, csfDest);

% IMPORTANT: load headers from the permanent copied files
GM = spm_vol(gmDest);
WM = spm_vol(wmDest);
CSF = spm_vol(csfDest);

% Now it is safe to delete the temporary SPM outputs
ea_delete(c1Path);
ea_delete(c2Path);
ea_delete(c3Path);

end

function [img, rng] = normalizeImage(img)
rng = prctile(img(:),[1,99,0,100]);
img = (img-rng(1))/(rng(2)-rng(1));
end

function [MRInfo, MRImg, CTInfo, CTImg] = loadReslicedMRCT(mrFile, ctFile)
% Reslice MR and CT to 0.4 x 0.4 x 0.4 mm using SPM & Lead-DBS helpers.
% Returns resliced volumes + headers, ready to assign to app.MRInfo/app.MRImg/etc.

% Desired isotropic resolution
newRes = [0.4 0.4 0.4];   % mm

% --- Load MR header
Vmr = spm_vol(mrFile);

% Use Lead-DBS helper to get voxel size robustly
niiMr    = ea_load_nii(mrFile);
oldResMr = double(niiMr.voxsize(:)');   % 1x3 [dx dy dz]

% Compute new dimensions so physical FOV is approximately preserved
scale  = oldResMr ./ newRes;           % how many more voxels along each axis
newDim = ceil(Vmr.dim .* scale);       % new [nx ny nz]

% --- Build a reference volume with desired geometry (0.4 mm)
[mrPath, mrName, mrExt] = fileparts(mrFile);
refPath = fullfile(mrPath, [mrName '_ref0p4' mrExt]);

Vref       = Vmr;
Vref.fname = refPath;
Vref.mat(1:3,1:3) = diag(newRes);  % set voxel size in affine
Vref.dim   = newDim;

% Create the empty reference NIfTI on disk
Vref = spm_create_vol(Vref);

% --- Set reslice options
flags = struct();
flags.mask   = 1;
flags.mean   = 0;
flags.interp = 1;       % trilinear (OK for MR/CT here)
flags.which  = 1;       % write resliced image only
flags.wrap   = [0 0 0];
flags.prefix = 'r_';

% ============================================================
% 1) Reslice MR into the 0.4 mm reference space
% ============================================================
Pmr = char(refPath, mrFile);   % first row = reference, second = moving
spm_reslice(Pmr, flags);

mrResliced = fullfile(mrPath, ['r_' mrName mrExt]);

% ============================================================
% 2) Reslice CT into the same 0.4 mm reference space
%    (assuming CT is already coregistered to MR)
% ============================================================
[ctPath, ctName, ctExt] = fileparts(ctFile);
Pct = char(refPath, ctFile);
spm_reslice(Pct, flags);

ctResliced = fullfile(ctPath, ['r_' ctName ctExt]);

% ============================================================
% 3) Load resliced volumes for use in the app
% ============================================================
MRInfo = spm_vol(mrResliced);
MRImg  = spm_read_vols(MRInfo);

CTInfo = spm_vol(ctResliced);
CTImg  = spm_read_vols(CTInfo);
end

function approved = isBrainshiftApproved_local(subjRoot, subjId)

% Returns true ONLY when the Lead-DBS brainshift method JSON exists
% and explicitly contains approval == 1.
%
% Missing brainshift folder / missing JSON / rejected / malformed JSON all return false.

approved = false;

jsonFile = fullfile( ...
    subjRoot, ...
    'brainshift', ...
    'log', ...
    [subjId '_desc-brainshiftmethod.json']);

% Brain shift was never run, or no log exists
if ~isfile(jsonFile)
    return;
end

try
    info = jsondecode(fileread(jsonFile));

    if isfield(info, 'approval') && ...
            isnumeric(info.approval) && ...
            isscalar(info.approval) && ...
            info.approval == 1

        approved = true;
    end

catch ME
    warning( ...
        'Could not read brain-shift approval file "%s": %s', ...
        jsonFile, ME.message);
end
end