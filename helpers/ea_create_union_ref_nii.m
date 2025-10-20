function ref = ea_create_union_ref_nii(outPath, inFiles, voxSize)
% CREATE_UNION_REF_NII  Build an empty reference NIfTI covering the union
% of multiple inputs, at a chosen voxel size.
%
% outPath : output .nii or .nii.gz for the reference
% inFiles : char or cellstr of input NIfTI paths
% voxSize : [vx vy vz] in mm (optional). If omitted/[], uses median of inputs.
%
% Returns struct 'ref' with fields:
%   .Mout   (4x4)  voxel->world affine used for the reference
%   .outSz  (1x3)  output image size
%   .vox    (1x3)  voxel size used
%   .bbox   (2x3)  [min; max] world-mm bounds (snapped to grid)
%   .path   (char) written file path (outPath)
%
% Requires SPM12 on MATLAB path.

if ischar(inFiles), inFiles = cellstr(inFiles); end
if nargin < 3, voxSize = []; end

% --- 1) Read input headers, collect world-space corners
N = numel(inFiles);
V = cell(N,1); M = cell(N,1); dim = zeros(N,3); vox = zeros(N,3);
allCorners = [];

for k = 1:N
    V{k}   = spm_vol(inFiles{k});
    M{k}   = V{k}.mat;                % voxel->world (mm)
    dim(k,:) = V{k}.dim(1:3);
    vox(k,:) = [norm(M{k}(1:3,1)) norm(M{k}(1:3,2)) norm(M{k}(1:3,3))];

    % 8 voxel-center corners
    corn = [
        1 1 1;
        dim(k,1) 1 1;
        1 dim(k,2) 1;
        1 1 dim(k,3);
        dim(k,1) dim(k,2) 1;
        dim(k,1) 1 dim(k,3);
        1 dim(k,2) dim(k,3);
        dim(k,1) dim(k,2) dim(k,3)];
    cw = (M{k} * [corn, ones(8,1)].').';
    allCorners = [allCorners; cw(:,1:3)]; %#ok<AGROW>
end

% --- 2) Decide voxel size and union bbox (snap to grid)
if isempty(voxSize)
    voxSize = median(vox,1);
end
vx = voxSize(1); vy = voxSize(2); vz = voxSize(3);

xmin = floor(min(allCorners(:,1))/vx)*vx; xmax = ceil(max(allCorners(:,1))/vx)*vx;
ymin = floor(min(allCorners(:,2))/vy)*vy; ymax = ceil(max(allCorners(:,2))/vy)*vy;
zmin = floor(min(allCorners(:,3))/vz)*vz; zmax = ceil(max(allCorners(:,3))/vz)*vz;

Xw = xmin:vx:xmax; Yw = ymin:vy:ymax; Zw = zmin:vz:zmax;
outSz = [numel(Xw) numel(Yw) numel(Zw)];

% --- 3) Build output affine (SPM: mm = M * [i j k 1]')
Mout = eye(4);
Mout(1,1) = vx; Mout(2,2) = vy; Mout(3,3) = vz;
Mout(1,4) = xmin - vx;   % i=1 -> xmin
Mout(2,4) = ymin - vy;   % j=1 -> ymin
Mout(3,4) = zmin - vz;   % k=1 -> zmin

% --- 4) Write empty reference (zeros) with SPM
wantGz = endsWith(outPath,'.gz','IgnoreCase',true);
niiPath = outPath;
if wantGz, niiPath = regexprep(outPath,'\.gz$','', 'ignorecase'); end
if ~endsWith(niiPath,'.nii','IgnoreCase',true), niiPath = [niiPath '.nii']; end

Vo = struct('fname', niiPath, ...
            'dim', outSz, ...
            'dt', [16 0], ...      % float32
            'mat', Mout, ...
            'pinfo', [1;0;0]);     % scale=1, offset=0
Vo = spm_create_vol(Vo);
spm_write_vol(Vo, zeros(outSz,'single'));   % zero-filled reference

if wantGz
    try
        gzip(niiPath);
        delete(niiPath);
    catch
        warning('gzip failed; leaving uncompressed NIfTI at %s', niiPath);
        wantGz = false;
    end
end

% --- 5) Return handy metadata
ref.Mout  = Mout;
ref.outSz = outSz;
ref.vox   = [vx vy vz];
ref.bbox  = [xmin ymin zmin; xmax ymax zmax];
ref.path  = outPath;
fprintf('Reference NIfTI written: %s\n', ref.path);
end