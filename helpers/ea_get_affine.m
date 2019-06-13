function best_affine = ea_get_affine(nii, type)
% Get the affine matrix of the NIfTI file
%
% Default type is 'SPM', the affine matrix is for one-based voxel to world
% space transformation. Otherwise, the affine matrix is used for zero-based
% calculation.

if nargin < 2
    type = 'SPM';
end

% Remove volume index used in SPM (',1' in '/PATH/TO/image.nii.gz,1')
[fpath, ~, fext] = ea_niifileparts(nii);

hdr = ea_fslhd([fpath, fext]);

if hdr.sform_code ~= 0 % Prefer sform
    affine = [hdr.sto_xyz1; hdr.sto_xyz2; hdr.sto_xyz3; hdr.sto_xyz4];
elseif hdr.qform_code ~= 0
    affine = [hdr.qto_xyz1; hdr.qto_xyz2; hdr.qto_xyz3; hdr.qto_xyz4];
else
    warning('Neither sform nor qform detected! Fallback to base affine.');
    ndims = hdr.dim0;
    shape = zeros(1, ndims);
    zooms = zeros(1, ndims);
    for i=1:ndims
        eval(['shape(i) = hdr.dim', num2str(i), ';']);
        eval(['zooms(i) = hdr.pixdim', num2str(i), ';']);
    end
    if ndims >= 3
        shape = shape(1:3);
        zooms = zooms(1:3);
    else
        full_shape = ones(1, 3);
        full_zooms = ones(1, 3);
        full_shape(1:ndims) = shape;
        full_zooms(1:ndims) = zooms;
        shape = full_shape;
        zooms = full_zooms;
    end
    zooms(1) = zooms(1) * -1; % Radiology orientation.
    % Get translations from center of image
    origin = (shape - 1) / 2.0;
    affine = eye(4);
    affine(1:3, 1:3) = diag(zooms);
    affine(1:3, end) = -origin .* zooms;
end

switch type
    case {'spm', 'SPM', 1, '1'}  % SPM type, one-based
        best_affine = affine;
        best_affine(:,4) = best_affine(:,4) - sum(best_affine(:,1:3),2);
    case {0, '0'}  % other types, zero-based
        best_affine = affine;
end
