function [mat, spmMat] = ea_ants_image_affine(fname)
% Voxel-to-RAS affine used by the bundled ANTs reader (SPM, one-based voxels).
% SPM prefers sform, whereas ITK may use qform, e.g. for a sheared CT sform.
% Ask ANTs which geometry it reads instead of duplicating ITK's version-
% dependent selection rules. Use its rounded text only to select the exact
% NIfTI affine, never as the returned geometry. No image/header is modified.

V = spm_vol(fname);
spmMat = V.mat;
sform = V.private.mat;
qform = readQform(V.private);
[x, y, z] = ndgrid([1 V.dim(1)], [1 V.dim(2)], [1 V.dim(3)]);
corners = [x(:) y(:) z(:) ones(8,1)]';
if ~isempty(V.private.mat0_intent) && ...
        max(abs((sform-qform)*corners), [], 'all') < 1e-4
    mat = sform;
    return;
end

ea_libs_helper;
header = ea_getExec(fullfile(ea_getearoot, 'ext_libs', 'ANTs', 'PrintHeader'));
[status, output] = ea_runcmd([header, ' ', ea_path_helper(fname)]);
spacing = regexp(output, '(?m)^\s*Spacing\s+\[([^\]]+)\]', 'tokens', 'once');
origin = regexp(output, '(?m)^\s*Origin\s+\[([^\]]+)\]', 'tokens', 'once');
direction = regexp(output, ...
    '(?m)^\s*Direction\s*\r?\n([^\r\n]+)\r?\n([^\r\n]+)\r?\n([^\r\n]+)', ...
    'tokens', 'once');
if status ~= 0 || isempty(spacing) || isempty(origin) || isempty(direction)
    error('LeadDBS:ANTsImageGeometryRead', ...
        'Could not read the ANTs image geometry for %s:\n%s', fname, output);
end
spacing = sscanf(strrep(spacing{1}, ',', ' '), '%f');
origin = sscanf(strrep(origin{1}, ',', ' '), '%f');
direction = sscanf(strjoin(direction, ' '), '%f');
if numel(spacing) ~= 3 || numel(origin) ~= 3 || numel(direction) ~= 9 || ...
        any(~isfinite([spacing; origin; direction])) || any(spacing <= 0)
    error('LeadDBS:ANTsImageGeometryRead', ...
        'Invalid ANTs image geometry for %s.', fname);
end
direction = reshape(direction, 3, 3)';

% PrintHeader reports LPS, with its origin at voxel [0 0 0]. Compare all
% three components separately at its six-significant-digit text precision.
candidates = {sform, qform};
scores = inf(1,2);
for k = 1:numel(candidates)
    candidate = candidates(k);
    affine = diag([-1 -1 1 1])*candidate{1}*[eye(3) ones(3,1); 0 0 0 1];
    candidateSpacing = vecnorm(affine(1:3,1:3))';
    candidateDirection = affine(1:3,1:3)./candidateSpacing';
    expected = [candidateSpacing; affine(1:3,4); candidateDirection(:)];
    observed = [spacing; origin; direction(:)];
    residual = abs(expected-observed)./max(1, abs(observed));
    if all(residual <= 6e-6)
        scores(k) = norm(residual);
    end
end
[score, best] = min(scores);
if isfinite(score)
    mat = candidates{best};
    return;
end
error('LeadDBS:ANTsImageGeometryMismatch', ...
    ['ANTs reads a geometry matching neither NIfTI sform nor qform for %s. ' ...
     'Cannot safely convert its CT registration to SPM coordinates.'], fname);
end

function mat = readQform(volume)
% Decode as nifti_quatern_to_mat44 in the NIfTI reference library used by
% ITK. SPM's mat0 uses a different near-180-degree threshold, so simply
% returning mat0 can introduce another rotation error for rounded headers.
% https://github.com/NIFTI-Imaging/nifti_clib/blob/master/niftilib/nifti1_io.c
info = struct(volume);
h = info.hdr;
if h.qform_code <= 0
    mat = volume.mat0;
    return;
end
q = double([h.quatern_b h.quatern_c h.quatern_d]);
aSquared = 1-sum(q.^2);
if aSquared < 1e-7
    q = q/norm(q);
    a = 0;
else
    a = sqrt(aSquared);
end
b = q(1); c = q(2); d = q(3);
rotation = [a*a+b*b-c*c-d*d, 2*(b*c-a*d), 2*(b*d+a*c); ...
            2*(b*c+a*d), a*a+c*c-b*b-d*d, 2*(c*d-a*b); ...
            2*(b*d-a*c), 2*(c*d+a*b), a*a+d*d-c*c-b*b];
spacing = double(h.pixdim(2:4));
spacing(spacing <= 0) = 1;
if h.pixdim(1) < 0
    spacing(3) = -spacing(3);
end
origin = double([h.qoffset_x; h.qoffset_y; h.qoffset_z]);
mat = double(single([rotation*diag(spacing), origin; 0 0 0 1]));
% Both returned affines use millimetres and one-based voxel coordinates.
switch bitand(h.xyzt_units,7)
    case 1 % metres
        mat(1:3,:) = 1000*mat(1:3,:);
    case 3 % micrometres
        mat(1:3,:) = 0.001*mat(1:3,:);
end
mat = mat*[eye(3) -ones(3,1); 0 0 0 1];
end
