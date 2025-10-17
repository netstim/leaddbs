function app = ea_prepare_legui(app, event) %#ok<INUSD>
    preferBrainshiftCT = true;
    enforceSameGrid    = false;
    useMNI             = false;
    app.LoadImgsBtnH.Enable = "off";

    % ------- Determine subject root (Lead-DBS subject folder) -------
    app.SelDir = LeG_lastDir();
    subjRoot = app.SelDir;
    [~, subjId] = fileparts(subjRoot);
    app.SelDir = subjRoot;
    
    GM_filename = strcat(subjId, '_ses-preop_space-anchorNative_desc-preproc_acq-iso_label-GM_mod-iso_T1w_mask.nii');
    WM_filename = strcat(subjId, '_ses-preop_space-anchorNative_desc-preproc_acq-iso_label-WM_mod-iso_T1w_mask.nii');
    Electrodes_filename = strcat(subjId, '_electrodes.mat');
    Channels_filename = strcat(subjId, '_channels.mat');
    % ------- Find MR/CT inside that Lead-DBS tree -------
    try
        [mrFile, ctFile] = findLeadDBSImagePair_local(subjRoot, preferBrainshiftCT, useMNI);
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

    [app.MRImg, MRRng] = normalizeImage(app.MRImg);
    [app.CTImg, app.CTRng] = normalizeImage(app.CTImg);
    
    if exist(fullfile(subjRoot, 'coregistration', 'anat', GM_filename), 'file') && ...
       exist(fullfile(subjRoot, 'coregistration', 'anat', WM_filename), 'file')
    
        % load existing GM/WM
        GM = spm_vol(fullfile(subjRoot, 'coregistration', 'anat', GM_filename));
        WM = spm_vol(fullfile(subjRoot, 'coregistration', 'anat', GM_filename));
    else
        % run segmentation
        [GM, WM] = get_segmentations(mrFile, subjRoot, GM_filename, WM_filename);
    end
    app.GrayInfo = GM;
    app.WhiteInfo = WM;
    app.GrayImg = spm_read_vols(app.GrayInfo);
    app.WhiteImg = spm_read_vols(app.WhiteInfo);

    app.BrainSurfRaw = []; app.ProjSurfRaw = [];
    [app.BrainSurfRaw,app.ProjSurfRaw] = LeG_genSurfaces({app.GrayImg,app.WhiteImg},app.GrayInfo);

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

    app.ElectFile = app.recofile;
    if isfile(app.ElectFile)
        ES = load(app.ElectFile);
        copyIf_local(app, ES, {'ElecXYZRaw','ElecFullIdxRaw','ElecCOMIdxRaw',...
                               'ElecXYZProjRaw','ElecFullIdxProjRaw','ElecCOMIdxProjRaw',...
                               'ElecMapRaw','DepthElecRaw','MicroElecRaw','RefElecRaw','GndElecRaw'});
        if isfield(ES,'PatientIDStr') && ~strcmp(app.PatientIDStr, ES.PatientIDStr)
            app.PatientIDStr    = ES.PatientIDStr;
            app.MainFigure.Name = sprintf('LeGUI v%0.1f - %s', app.Version, app.PatientIDStr);
        end
    end
end

% -------------------- helpers (local scope) --------------------

function [mrFile, ctFile] = findLeadDBSImagePair_local(rootDir, preferBrainshift, useMNI)
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
        if preferBrainshift
            ctPats = [ ...
              "brainshift/*_ses-postop_space-anchorNative_rec-tonemappedbrainshift_desc-preproc_CT.nii", ...
              "coregistration/anat*/**/*_ses-postop_space-anchorNative_*desc-preproc_CT.nii", ...
              "preprocessing/*_ses-postop_space-anchorNative_desc-preproc_CT.nii" ...
            ];
        else
            ctPats = [ ...
              "coregistration/anat*/**/*_ses-postop_space-anchorNative_*desc-preproc_CT.nii", ...
              "preprocessing/*_ses-postop_space-anchorNative_desc-preproc_CT.nii", ...
              "brainshift/*_ses-postop_space-anchorNative_rec-tonemappedbrainshift_desc-preproc_CT.nii" ...
            ];
        end
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

function [GM, WM] = get_segmentations(mrPath, subjRoot, GM_filename, WM_filename)
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

    GM = spm_vol(c1Path); %c1Vol = spm_read_vols(V1);
    WM = spm_vol(c2Path); %c2Vol = spm_read_vols(V2);
    copyfile(c1Path, fullfile(subjRoot, 'coregistration', 'anat', GM_filename));
    copyfile(c2Path, fullfile(subjRoot, 'coregistration', 'anat', WM_filename));
end

function [img, rng] = normalizeImage(img)
    rng = prctile(img(:),[1,99,0,100]);
    img = (img-rng(1))/(rng(2)-rng(1));
end