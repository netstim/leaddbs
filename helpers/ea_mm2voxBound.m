function [vox, removed] = ea_mm2voxBound(mm, reference, dim, base)
% Converts mm-coordinates to voxel-coordinates and bound the coordinates
% within the image dimension.
% coords need to be row vector: N*3

if ischar(reference)
    if ~exist('base', 'var')
        base = 1; % Use one-based indexing by default
    end
    affine = ea_get_affine(reference, base);
elseif isnumeric(reference)
    if ~exist('dim', 'var') || isempty(dim)
        ea_cprintf('CmdWinWarnings', 'Reference image dimension has to be defined to bound the voxel coordinates!\n');
    end
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
    removed = ~prod(vox>=1,2);
else % zero-based indexing
    removed = ~prod(vox>=0,2);
end

% Filter voxel coordinates outside of the image dimension
if exist('dim', 'var') && ~isempty(dim)
    removed = removed | vox(:,1)>dim(1) | vox(:,2)>dim(2) | vox(:,3)>dim(3);
elseif ischar(reference)
    header = ea_fslhd(reference);
    removed = removed | vox(:,1)>header.dim1 | vox(:,2)>header.dim2 | vox(:,3)>header.dim3;
end

vox(removed, :) = [];

