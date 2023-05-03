function vox = ea_mm2voxBound(mm, reference, base)
% Converts mm-coordinates to voxel-coordinates and bound the coordinates
% within the image dimension.
% coords need to be row vector: N*3

if ischar(reference)
    if ~exist('base', 'var')
        base = 1; % Use one-based indexing by default
    end
    affine = ea_get_affine(reference, base);
elseif isnumeric(reference)
    ea_cprintf('CmdWinWarnings', 'Reference image has to be defined to bound the voxel coordinates!\n');
    affine = reference;
else
    ea_cprintf('CmdWinErrors', 'Reference image has to be defined to bound the voxel coordinates!\n');
    return;
end

vox = [mm, ones(size(mm,1),1)] / affine';
vox(:,4) = [];

vox = round(vox);

% Filter zero/negative voxel coordinates
if ~exist('base', 'var') || base % one-based indexing
    vox(~prod(vox>=1,2), :) = [];
else % zero-based indexing
    vox(~prod(vox>=0,2), :) = [];
end

% Filter voxel coordinates outside of the image dimension
if ischar(reference)
    header = ea_fslhd(reference);
    vox(vox(:,1)>header.dim1 | vox(:,2)>header.dim2 | vox(:,3)>header.dim3, :) = [];
end
