function ea_diagnose_structural_CM_data(subjectDir, parcellationName)
% EA_DIAGNOSE_STRUCTURAL_CM_DATA  Check data that feeds into structural connectivity matrix.
% Reads and reports only - does not modify any Lead-DBS code or pipeline.
%
% Usage:
%   ea_diagnose_structural_CM_data(subjectDir, parcellationName)
%
% Example:
%   ea_diagnose_structural_CM_data(...
%     '/Users/maikemustin/Desktop/BIDSTest/derivatives/leaddbs/sub-patienttestlegacy', ...
%     'AAL')
%
% If parcellationName is omitted, the script will try to find a b0w*.nii in
% templates/labeling/ and use that name.

if nargin < 1 || isempty(subjectDir)
    subjectDir = uigetdir([], 'Select subject directory (e.g. .../derivatives/leaddbs/sub-XXX)');
    if isequal(subjectDir, 0), return; end
end
if subjectDir(end) == filesep, subjectDir = subjectDir(1:end-1); end

if nargin < 2 || isempty(parcellationName)
    labDir = fullfile(subjectDir, 'templates', 'labeling');
    b0wFiles = dir(fullfile(labDir, 'b0w*.nii'));
    if isempty(b0wFiles)
        fprintf('No b0w* parcellation found in %s\n', labDir);
        parcellationName = '';
    else
        [~, name] = fileparts(b0wFiles(1).name);
        parcellationName = strrep(name, 'b0w', '');
        fprintf('Auto-detected parcellation from file: %s\n', b0wFiles(1).name);
    end
end
% If given name not found, try to find any b0w*.nii and use it
labDir = fullfile(subjectDir, 'templates', 'labeling');
b0wFiles = dir(fullfile(labDir, 'b0w*.nii'));
if ~isempty(parcellationName) && ~isfile(fullfile(labDir, ['b0w', parcellationName, '.nii'])) && ~isempty(b0wFiles)
    [~, name] = fileparts(b0wFiles(1).name);
    foundName = strrep(name, 'b0w', '');
    fprintf('Note: b0w%s.nii not found; using instead: b0w%s.nii\n', parcellationName, foundName);
    parcellationName = foundName;
end

fprintf('\n=== Structural CM data diagnostic (read-only) ===\n');
fprintf('Subject dir: %s\n', subjectDir);
fprintf('Parcellation: %s\n\n', parcellationName);

%% 1) B0 image
b0Path = fullfile(subjectDir, 'preprocessing', 'dwi');
d = dir(fullfile(b0Path, '*_b0.nii'));
if isempty(d), d = dir(fullfile(b0Path, '*b0*.nii')); end
if isempty(d)
    fprintf('[1] B0: NOT FOUND in %s\n', b0Path);
    Vb0 = [];
else
    b0File = fullfile(d(1).folder, d(1).name);
    Vb0 = spm_vol(b0File);
    fprintf('[1] B0: %s\n', d(1).name);
    fprintf('    dim = [%d %d %d]\n', Vb0.dim);
    fprintf('    mat(1:3,4) (origin approx) = [%.2f %.2f %.2f]\n', Vb0.mat(1:3,4));
end

%% 2) Parcellation in b0 space (b0w*parcellation*.nii)
labDir = fullfile(subjectDir, 'templates', 'labeling');
parcFile = fullfile(labDir, ['b0w', parcellationName, '.nii']);
if ~isfile(parcFile)
    fprintf('[2] Parcellation (b0w): NOT FOUND: %s\n', parcFile);
    Vatl = [];
else
    Vatl = spm_vol(parcFile);
    parcImg = spm_read_vols(Vatl);
    nz = parcImg(parcImg(:) ~= 0);
    fprintf('[2] Parcellation (b0w): %s\n', ['b0w', parcellationName, '.nii']);
    fprintf('    dim = [%d %d %d]\n', Vatl.dim);
    fprintf('    non-zero voxels: %d (total %d)\n', numel(nz), numel(parcImg));
    if isempty(nz)
        fprintf('    *** WARNING: Parcellation is ALL ZERO -> CM will be zero.\n');
    else
        u = unique(nz);
        fprintf('    unique labels (first 20): %s\n', mat2str(u(1:min(20,numel(u)))'));
    end
    if ~isempty(Vb0) && (any(Vatl.dim ~= Vb0.dim) || norm(Vatl.mat - Vb0.mat) > 1e-6)
        fprintf('    *** WARNING: dim or mat differ from B0 -> coordinate mismatch possible.\n');
    end
end

%% 3) FTR (fiber tract)
ftrPath = fullfile(subjectDir, 'connectomics', 'dMRI', 'FTR.mat');
if ~isfile(ftrPath)
    fprintf('[3] FTR: NOT FOUND: %s\n', ftrPath);
    fibs = []; idx = []; voxmm = ''; ftrMat = [];
else
    f = load(ftrPath);
    fibs = f.fibers;
    idx = f.idx;
    if isfield(f, 'voxmm'), voxmm = f.voxmm; else, voxmm = '?'; end
    if isfield(f, 'mat'), ftrMat = f.mat; else, ftrMat = []; end
    fprintf('[3] FTR: connectomics/dMRI/FTR.mat\n');
    fprintf('    fibers size: [%d x %d]\n', size(fibs,1), size(fibs,2));
    fprintf('    voxmm: %s\n', voxmm);
    fprintf('    fiber coord range: X [%.2f, %.2f], Y [%.2f, %.2f], Z [%.2f, %.2f]\n', ...
        min(fibs(:,1)), max(fibs(:,1)), min(fibs(:,2)), max(fibs(:,2)), min(fibs(:,3)), max(fibs(:,3)));
    if ~isempty(Vb0)
        % Expected voxel range for b0 (1-based)
        fprintf('    B0 voxel range: X [1, %d], Y [1, %d], Z [1, %d]\n', Vb0.dim(1), Vb0.dim(2), Vb0.dim(3));
        if strcmp(voxmm, 'vox')
            if max(fibs(:,1)) > Vb0.dim(1) || max(fibs(:,2)) > Vb0.dim(2) || max(fibs(:,3)) > Vb0.dim(3)
                fprintf('    *** WARNING: Fiber voxel coords exceed B0 dims.\n');
            end
        end
    end
end

%% 4) What would spm_sample_vol return?
% SPM's spm_sample_vol expects coordinates in MM (world space), not voxel.
% If FTR stores voxel and the pipeline passes them without conversion,
% sampling would use wrong positions (voxel numbers interpreted as mm).
if ~isempty(Vatl) && ~isempty(fibs) && size(fibs,1) > 0
    fprintf('\n[4] Sampling parcellation at fiber coordinates:\n');
    if strcmpi(voxmm, 'vox')
        % Direkte Voxel-Sampling-Variante (entspricht aktueller CM-Berechnung)
        nSample = min(1000, size(fibs,1));
        parcImg = spm_read_vols(Vatl);
        x = round(fibs(1:nSample,1));
        y = round(fibs(1:nSample,2));
        z = round(fibs(1:nSample,3));
        dim = Vatl.dim;
        inBounds = x >= 1 & x <= dim(1) & ...
                   y >= 1 & y <= dim(2) & ...
                   z >= 1 & z <= dim(3);
        linIdx = sub2ind(dim, x(inBounds), y(inBounds), z(inBounds));
        vals = zeros(nSample,1);
        vals(inBounds) = parcImg(linIdx);
        nHit = sum(vals(:) ~= 0 & ~isnan(vals(:)));
        fprintf('    (Assuming fibers are VOXEL; sampled via direct voxel indexing.)\n');
        fprintf('    Sample of first %d points: %d non-zero parcellation labels.\n', min(1000, size(fibs,1)), nHit);
        if nHit == 0
            fprintf('    *** If pipeline passes voxel coords to spm_sample_vol (which expects mm),\n');
            fprintf('       all samples would be at wrong positions -> CM zero.\n');
        end
    else
        % Assume mm
        nS = min(1000, size(fibs,1));
        sample = spm_sample_vol(Vatl, double(fibs(1:nS,1)), double(fibs(1:nS,2)), double(fibs(1:nS,3)), 0);
        nHit = sum(sample(:) ~= 0 & ~isnan(sample(:)));
        fprintf('    (Assuming fibers are MM.)\n');
        fprintf('    Sample of first %d points: %d non-zero parcellation labels.\n', min(1000, size(fibs,1)), nHit);
    end
end

%% 5) Anatomical image (preprocessed T1)
anatDir = fullfile(subjectDir, 'preprocessing', 'anat');
Va = [];
if ~isfolder(anatDir)
    fprintf('\n[5] Anatomical (preproc T1): preprocessing/anat directory NOT FOUND.\n');
else
    % Try to find a reasonably named T1 (BIDS-style)
    cand = dir(fullfile(anatDir, '*desc-preproc*_T1w.nii'));
    if isempty(cand)
        cand = dir(fullfile(anatDir, '*acq-iso_T1w.nii'));
    end
    if isempty(cand)
        cand = dir(fullfile(anatDir, '*T1w.nii'));
    end
    if isempty(cand)
        fprintf('\n[5] Anatomical (preproc T1): NOT FOUND in %s\n', anatDir);
    else
        aFile = fullfile(cand(1).folder, cand(1).name);
        Va = spm_vol(aFile);
        fprintf('\n[5] Anatomical (preproc T1): %s\n', cand(1).name);
        fprintf('    dim = [%d %d %d]\n', Va.dim);
        fprintf('    mat(1:3,4) (origin approx) = [%.2f %.2f %.2f]\n', Va.mat(1:3,4));
        if ~isempty(Vb0)
            dMat = Va.mat - Vb0.mat;
            fprintf('    |T1.mat - B0.mat|_F = %.4f\n', norm(dMat(:)));
        end
    end
end

%% 6) Normalization file (y_ea_inv_normparams.nii)
normFile = fullfile(subjectDir, 'y_ea_inv_normparams.nii');
Vn = [];
if isfile(normFile)
    try
        Vn = spm_vol(normFile);
        fprintf('\n[6] Normalization file: y_ea_inv_normparams.nii\n');
        fprintf('    dim = [%d %d %d]\n', Vn.dim);
        fprintf('    mat(1:3,4) = [%.2f %.2f %.2f]\n', Vn.mat(1:3,4));
    catch ME
        fprintf('\n[6] Normalization file: error reading %s\n    %s\n', normFile, ME.message);
    end
else
    fprintf('\n[6] Normalization file: NOT FOUND: %s\n', normFile);
end

%% 7) Coregistration transforms between B0 and T1
fprintf('\n[7] Coregistration transforms (B0 <-> T1):\n');
coregDirs = {
    fullfile(subjectDir, 'coregistration', 'transformations')
    fullfile(subjectDir, 'preprocessing', 'anat')
    subjectDir
};
if isempty(Vb0) || isempty(Va)
    fprintf('    (Skipping detailed search because B0 or T1 was not found.)\n');
else
    [~, b0Name] = fileparts(Vb0.fname);
    [~, anatName] = fileparts(Va.fname);
    b0Name = regexprep(b0Name, '\.nii(\.gz)?$', '');
    anatName = regexprep(anatName, '\.nii(\.gz)?$', '');
    pattern1 = [b0Name, '2', anatName]; % b0 -> T1
    pattern2 = [anatName, '2', b0Name]; % T1 -> b0

    foundAny = false;
    for iDir = 1:numel(coregDirs)
        cdir = coregDirs{iDir};
        if ~isfolder(cdir), continue; end
        mfiles = dir(fullfile(cdir, '*.mat'));
        if isempty(mfiles), continue; end
        for k = 1:numel(mfiles)
            fn = mfiles(k).name;
            if contains(fn, pattern1) || contains(fn, pattern2)
                foundAny = true;
                fprintf('    %s (in %s)\n', fn, cdir);
            end
        end
    end
    if ~foundAny
        fprintf('    No explicit *.mat transforms with patterns "%s" or "%s" found.\n', pattern1, pattern2);
    end
end

%% 8) FTR affine vs. B0 / Parcellation
if ~isempty(fibs) && ~isempty(Vb0)
    fprintf('\n[8] FTR affine vs. B0 / Parcellation:\n');
    if ~isempty(ftrMat)
        fprintf('    FTR.mat present. ||FTR.mat - B0.mat||_F = %.4f\n', norm(ftrMat(:) - Vb0.mat(:)));
        if ~isempty(Vatl)
            fprintf('    ||FTR.mat - b0wParc.mat||_F = %.4f\n', norm(ftrMat(:) - Vatl.mat(:)));
        end
    else
        fprintf('    FTR.mat does not contain an affine "mat" field.\n');
    end
end

fprintf('\n=== End diagnostic ===\n\n');
