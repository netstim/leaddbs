function ref = ea_create_threshold_union_ref_nii(outPath, inFiles, thr, voxSize)
% CREATE_THRESHOLD_UNION_REF_NII
% Make a zero-filled reference NIfTI whose FOV is the union bounding box
% of all *thresholded* regions (V > thr) across input images.
%
% outPath : '.nii' or '.nii.gz'
% inFiles : char or cellstr of NIfTI paths
% thr     : numeric threshold (apply as V > thr in each native image)
% voxSize : [vx vy vz] mm (optional). If [], uses median voxel size across inputs.
%
% Returns: struct with .Mout .outSz .vox .bbox .path
%
% Requires: Lead-DBS (ea_load_nii / ea_write_nii) or SPM12 (spm_*).
% This version uses SPM to write; easy to switch to ea_write_nii (see comment).

if ischar(inFiles), inFiles = cellstr(inFiles); end
if nargin < 3 || isempty(thr), error('Provide a numeric threshold "thr".'); end
if nargin < 4, voxSize = []; end

N = numel(inFiles);
M = cell(N,1);
dim = zeros(N,3);
vox = zeros(N,3);
allCorners = [];

% --- scan each image just to find thresholded extents ---
for k = 1:N
    nii = ea_load_nii(inFiles{k});  % .img (double/single), .mat (voxel->world)
    Vk  = nii.img;
    M{k}= nii.mat;
    dim(k,:) = size(Vk);
    vox(k,:) = [norm(M{k}(1:3,1)) norm(M{k}(1:3,2)) norm(M{k}(1:3,3))];

    % threshold (avoid NaNs)
    mask = Vk > thr;
    mask(isnan(Vk)) = false;

    if ~any(mask(:))
        % no above-threshold voxels in this file -> skip
        continue
    end

    % Find min/max indices along each axis without full find()
    i_any = squeeze(any(any(mask,2),3));
    j_any = squeeze(any(any(mask,1),3));
    k_any = squeeze(any(any(mask,1),2));
    imin = find(i_any,1,'first'); imax = find(i_any,1,'last');
    jmin = find(j_any,1,'first'); jmax = find(j_any,1,'last');
    kmin = find(k_any,1,'first'); kmax = find(k_any,1,'last');

    % 8 voxel-center corners of the thresholded subvolume
    corn = [
        imin jmin kmin;
        imax jmin kmin;
        imin jmax kmin;
        imin jmin kmax;
        imax jmax kmin;
        imax jmin kmax;
        imin jmax kmax;
        imax jmax kmax];
    cw = (M{k} * [corn, ones(8,1)].').';   % world-mm corners (8x4)
    allCorners = [allCorners; cw(:,1:3)]; %#ok<AGROW>
end

if isempty(allCorners)
    error('No voxels exceeded the threshold in any input.');
end

% --- choose voxel size & union bbox (snap to grid) ---
if isempty(voxSize), voxSize = median(vox,1); end
vx = voxSize(1); vy = voxSize(2); vz = voxSize(3);

xmin = floor(min(allCorners(:,1))/vx)*vx; xmax = ceil(max(allCorners(:,1))/vx)*vx;
ymin = floor(min(allCorners(:,2))/vy)*vy; ymax = ceil(max(allCorners(:,2))/vy)*vy;
zmin = floor(min(allCorners(:,3))/vz)*vz; zmax = ceil(max(allCorners(:,3))/vz)*vz;

Xw = xmin:vx:xmax; Yw = ymin:vy:ymax; Zw = zmin:vz:zmax;
outSz = [numel(Xw) numel(Yw) numel(Zw)];

% --- build affine (SPM/Lead-DBS: mm = M * [i j k 1]') ---
Mout = eye(4);
Mout(1,1)=vx; Mout(2,2)=vy; Mout(3,3)=vz;
Mout(1,4)=xmin - vx; Mout(2,4)=ymin - vy; Mout(3,4)=zmin - vz;

% --- write zero-filled reference using SPM (fast for later spm_reslice) ---
wantGz = endsWith(outPath,'.gz','IgnoreCase',true);
niiPath = outPath;
if wantGz, niiPath = regexprep(outPath,'\.gz$','', 'ignorecase'); end
if ~endsWith(niiPath,'.nii','IgnoreCase',true), niiPath = [niiPath '.nii']; end

Vo = struct('fname', niiPath, 'dim', outSz, 'dt', [16 0], 'mat', Mout, 'pinfo', [1;0;0]);
Vo = spm_create_vol(Vo);
spm_write_vol(Vo, zeros(outSz,'single'));

if wantGz
    try gzip(niiPath); delete(niiPath); catch, warning('gzip failed; left %s', niiPath); end
end

% If you prefer Lead-DBS I/O, swap the write block with:
% tmpl = ea_load_nii(inFiles{1}); tmpl.img = zeros(outSz,'single'); tmpl.mat = Mout;
% tmpl.hdr.dime.dim(2:4)=outSz; tmpl.hdr.dime.pixdim(2:4)=[vx vy vz];
% tmpl.hdr.dime.datatype=16; tmpl.hdr.dime.bitpix=32;
% tmpl.hdr.hist.sform_code=1; tmpl.hdr.hist.qform_code=0;
% tmpl.hdr.hist.srow_x=Mout(1,:); tmpl.hdr.hist.srow_y=Mout(2,:); tmpl.hdr.hist.srow_z=Mout(3,:);
% ea_write_nii(tmpl, outPath);

% --- return metadata ---
ref.Mout  = Mout;
ref.outSz = outSz;
ref.vox   = [vx vy vz];
ref.bbox  = [xmin ymin zmin; xmax ymax zmax];
ref.path  = outPath;
fprintf('Thresholded reference written: %s\n', ref.path);
end