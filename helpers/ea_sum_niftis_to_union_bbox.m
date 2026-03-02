function outPath = ea_sum_niftis_to_union_bbox(outPath, inFiles, voxSize, interpMode)
% Sum multiple NIfTIs into a union-bounding-box reference (Lead-DBS friendly).
% Uses each file's .mat to respect orientation, different FOVs, etc.
%
% outPath    : 'sum.nii' or 'sum.nii.gz'
% inFiles    : char or cellstr of NIfTI paths
% voxSize    : [vx vy vz] mm (default: median of inputs)
% interpMode : 'linear' (default) or 'nearest' (use for masks)

if ischar(inFiles), inFiles = cellstr(inFiles); end
if nargin < 3 || isempty(voxSize), voxSize = []; end
if nargin < 4 || isempty(interpMode), interpMode = 'linear'; end
assert(ismember(interpMode, {'linear','nearest'}), 'interpMode must be linear|nearest');

% ---------- 1) Load headers, get corners ----------
N = numel(inFiles);
nii = cell(N,1); M = cell(N,1); dim = zeros(N,3); vox = zeros(N,3);
allCorners = [];

for k = 1:N
    nii{k} = ea_load_nii(inFiles{k});  % .img, .mat, .hdr
    M{k}   = nii{k}.mat;
    dim(k,:) = size(nii{k}.img);
    vox(k,:) = [norm(M{k}(1:3,1)) norm(M{k}(1:3,2)) norm(M{k}(1:3,3))];
    % 8 corners in voxel indices (centers)
    corn = [ ...
        1 1 1; dim(k,1) 1 1; 1 dim(k,2) 1; 1 1 dim(k,3); ...
        dim(k,1) dim(k,2) 1; dim(k,1) 1 dim(k,3); 1 dim(k,2) dim(k,3); ...
        dim(k,1) dim(k,2) dim(k,3)];
    cw = (M{k} * [corn, ones(8,1)].').';
    allCorners = [allCorners; cw(:,1:3)]; %#ok<AGROW>
end

% ---------- 2) Union bbox + output grid ----------
if isempty(voxSize), voxSize = median(vox,1); end
vx = voxSize(1); vy = voxSize(2); vz = voxSize(3);

xmin = floor(min(allCorners(:,1))/vx)*vx; xmax = ceil(max(allCorners(:,1))/vx)*vx;
ymin = floor(min(allCorners(:,2))/vy)*vy; ymax = ceil(max(allCorners(:,2))/vy)*vy;
zmin = floor(min(allCorners(:,3))/vz)*vz; zmax = ceil(max(allCorners(:,3))/vz)*vz;

Xw = xmin:vx:xmax; Yw = ymin:vy:ymax; Zw = zmin:vz:zmax;
[XX,YY,ZZ] = ndgrid(Xw, Yw, Zw);
outSz = [numel(Xw) numel(Yw) numel(Zw)];

% ---------- 3) Accumulate SUM ----------
sumVol = zeros(outSz, 'single');   % single for speed/bandwidth
invM = cellfun(@inv, M, 'uni', 0); % precompute inverses

for k = 1:N
    Vk = single(nii{k}.img);
    dk = dim(k,:);
    % ---- compute xi,yj,zk (world -> voxel indices for this Vk) ----
    nvox = numel(XX);
    pts  = [single(XX(:)) single(YY(:)) single(ZZ(:)) ones(nvox,1,'single')] * single(invM{k}).';
    xi   = reshape(pts(:,1), outSz);     % matches ndgrid shape of XX,YY,ZZ
    yj   = reshape(pts(:,2), outSz);
    zk   = reshape(pts(:,3), outSz);

    % ---- fast interpolant from grid VECTORS (avoids NDGRID-structure error) ----
    F = griddedInterpolant( ...
            {single(1:dk(1)), single(1:dk(2)), single(1:dk(3))}, ...
            Vk, interpMode, 'none');   % 'none' => NaN outside
    Vk_res = F(xi, yj, zk);
    Vk_res(isnan(Vk_res)) = 0;         % outside FOV -> 0

    sumVol = sumVol + Vk_res;
end

% ---------- 4) Output affine (SPM/Lead-DBS: mm = M * [i j k 1]') ----------
Mout = eye(4,'single');
Mout(1,1) = vx; Mout(2,2) = vy; Mout(3,3) = vz;
Mout(1,4) = xmin - vx; Mout(2,4) = ymin - vy; Mout(3,4) = zmin - vz;

% ---------- 5) Write with Lead-DBS (ea_write_nii) ----------
sumVol = single(sumVol);                        % class must match header (float32)

outnii           = nii{1};                      % template from first input
outnii.img       = sumVol;
outnii.mat       = double(Mout);

% -- keep header and .dim consistent with img --
outSz            = size(sumVol);                % [X Y Z]
outnii.dim       = outSz;                       % some Lead-DBS funcs read this
outnii.hdr.dime.dim(1)   = 3;                   % 3D image
outnii.hdr.dime.dim(2:4) = outSz;               % MUST equal size(img)
outnii.hdr.dime.dim(5:8) = 1;                   % clear extra dims
outnii.hdr.dime.pixdim(1)   = 1;                % qfac (leave 1)
outnii.hdr.dime.pixdim(2:4) = [vx vy vz];
outnii.hdr.dime.datatype = 16;                  % float32
outnii.hdr.dime.bitpix   = 32;

% -- sform (preferred) --
outnii.hdr.hist.sform_code = 1;
outnii.hdr.hist.qform_code = 0;
outnii.hdr.hist.srow_x = outnii.mat(1,:);
outnii.hdr.hist.srow_y = outnii.mat(2,:);
outnii.hdr.hist.srow_z = outnii.mat(3,:);

% IMPORTANT: call with the filename
outnii.fname=outPath;
ea_write_nii(outnii);
fprintf('Wrote summed NIfTI to %s\n', outPath);