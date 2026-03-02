function ea_visualize_structural_alignment(subjectDir)
% EA_VISUALIZE_STRUCTURAL_ALIGNMENT
% Visualisiert b0, b0w-Parzellierung, T1 und Fasern gemeinsam, um
% Ausrichtungs- und Affine-Probleme sichtbar zu machen.
%
% Aufruf (im MATLAB Command Window):
%   ea_visualize_structural_alignment( ...
%     '/Users/maikemustin/Desktop/BIDIBITS/derivatives/leaddbs/sub-patienttestlegacy');

if nargin < 1 || isempty(subjectDir)
    subjectDir = uigetdir([], 'Select subject directory (..../derivatives/leaddbs/sub-XXX)');
    if isequal(subjectDir,0), return; end
end
if subjectDir(end) == filesep, subjectDir = subjectDir(1:end-1); end

%% 1. Dateien finden
fprintf('\n=== Visual alignment diagnostic ===\n');
fprintf('Subject dir: %s\n', subjectDir);

% B0
b0Dir  = fullfile(subjectDir, 'preprocessing', 'dwi');
b0Cand = dir(fullfile(b0Dir, '*_b0.nii'));
if isempty(b0Cand)
    b0Cand = dir(fullfile(b0Dir, '*b0*.nii'));
end
if isempty(b0Cand)
    error('No B0 image found in %s', b0Dir);
end
b0File = fullfile(b0Cand(1).folder, b0Cand(1).name);
fprintf('B0: %s\n', b0Cand(1).name);

% Parzellierung im b0-Space (b0w*)
labDir = fullfile(subjectDir, 'templates', 'labeling');
parcCand = dir(fullfile(labDir, 'b0w*.nii'));
if isempty(parcCand)
    error('No b0w* parcellation found in %s', labDir);
end
parcFile = fullfile(parcCand(1).folder, parcCand(1).name);
fprintf('Parcellation (b0w): %s\n', parcCand(1).name);

% T1 (preproc)
anatDir = fullfile(subjectDir, 'preprocessing', 'anat');
t1Cand  = dir(fullfile(anatDir, '*desc-preproc*_T1w.nii'));
if isempty(t1Cand)
    t1Cand = dir(fullfile(anatDir, '*acq-iso_T1w.nii'));
end
if isempty(t1Cand)
    t1Cand = dir(fullfile(anatDir, '*T1w.nii'));
end
if isempty(t1Cand)
    error('No preproc T1 found in %s', anatDir);
end
t1File = fullfile(t1Cand(1).folder, t1Cand(1).name);
fprintf('T1 (preproc): %s\n', t1Cand(1).name);

% FTR
ftrFile = fullfile(subjectDir, 'connectomics', 'dMRI', 'FTR.mat');
if ~isfile(ftrFile)
    error('FTR.mat not found at %s', ftrFile);
end
fprintf('FTR: FTR.mat\n');

%% 2. NIfTIs laden
Vb0   = spm_vol(b0File);
b0img = spm_read_vols(Vb0);

Vparc   = spm_vol(parcFile);
parcImg = spm_read_vols(Vparc);

Vt1   = spm_vol(t1File);
t1img = spm_read_vols(Vt1);

ftr   = load(ftrFile);
fibs  = ftr.fibers;   % [N x 3], meist mm
idx   = ftr.idx;

fprintf('B0 dim = [%d %d %d], mat(1:3,4) = [%.2f %.2f %.2f]\n', ...
    Vb0.dim, Vb0.mat(1:3,4));
fprintf('b0w dim = [%d %d %d], non-zero voxels = %d\n', ...
    Vparc.dim, nnz(parcImg));

%% 3. b0 + b0wAtlas: 2D-Overlay (axial)
sliceZ = round(size(b0img,3)/2);  % mittlerer Slice
b0sl   = squeeze(b0img(:,:,sliceZ));
parcsl = squeeze(parcImg(:,:,sliceZ));
b0sl   = b0sl ./ (max(b0sl(:)) + eps);
parcMask = parcsl ~= 0;

figure('Name','b0 + b0wAtlas (axial)','Color','w');
imagesc(b0sl); colormap(gray); axis image off; hold on;
h = imagesc(parcMask);
set(h, 'AlphaData', 0.4);
colormap(gca, [gray(256); 1 1 0]);   % grau + gelb
title(sprintf('b0 (grau) + b0wAtlas (gelb), Slice Z=%d', sliceZ));

%% 4. T1 (auf b0-Gitter) + b0wAtlas
% Reslice T1 nach b0-Gitter
spm_reslice(char(Vb0.fname, Vt1.fname), struct('interp',1,'which',1,'mean',0));
[~,t1name,t1ext] = fileparts(Vt1.fname);
rt1File = fullfile(fileparts(Vt1.fname), ['r', t1name, t1ext]);
Vt1b    = spm_vol(rt1File);
t1b     = spm_read_vols(Vt1b);

if any(size(t1b) ~= size(b0img))
    warning('Resliced T1 size does not match B0 size, skipping T1 overlay.');
else
    t1sl = squeeze(t1b(:,:,sliceZ));
    figure('Name','T1 (auf b0-Gitter) + b0wAtlas','Color','w');
    imagesc(t1sl); colormap(gray); axis image off; hold on;
    h2 = imagesc(parcMask);
    set(h2, 'AlphaData', 0.4);
    colormap(gca, [gray(256); 1 1 0]);
    title(sprintf('T1 (grau) + b0wAtlas (gelb), Slice Z=%d', sliceZ));
end

%% 5. Fasern im b0-Voxelraum vs. b0 & Atlas (3D)
% Aktuell liegen die Fasern in FTR.mat bereits im VOXEL-Raum des b0
% (voxmm = 'vox', Affine = B0.mat). Daher hier KEINE weitere
% Transformation, sondern direkte Verwendung der Voxelkoordinaten.
XYZ_vox = fibs(:,1:3);

% B0- und Atlas-Maskenpunkte
b0mask   = b0img > max(b0img(:))/4;
[bx,by,bz] = ind2sub(size(b0img), find(b0mask));
parcmask = parcImg ~= 0;
[px,py,pz] = ind2sub(size(parcImg), find(parcmask));

figure('Name','3D: b0 / Atlas / Fasern','Color','w'); hold on;

% b0-FOV (hellgrau, stark ausgedünnt)
plot3(bx(1:500:end), by(1:500:end), bz(1:500:end), '.', ...
      'Color', [0.85 0.85 0.85], 'MarkerSize', 4);

% Atlas-Voxel (gelb)
plot3(px(1:50:end), py(1:50:end), pz(1:50:end), '.', ...
      'Color', [1 0.9 0], 'MarkerSize', 6);

% Fasern (rot)
nFibPlot = min(5000, size(XYZ_vox,1));
plot3(XYZ_vox(1:nFibPlot,1), XYZ_vox(1:nFibPlot,2), XYZ_vox(1:nFibPlot,3), ...
      'r.', 'MarkerSize', 4);

axis equal vis3d; view(3);
xlabel X; ylabel Y; zlabel Z;
title('b0-FOV (grau), b0wAtlas (gelb), Fasern (rot)');
legend({'b0','Atlas','Fasern'});

fprintf('=== End visual alignment diagnostic ===\n\n');
end
