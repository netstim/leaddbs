function [vox_overlap, mm_overlap, normVTAOverlap, normAtlasOverlap, mm_vta, mm_atlas]  = ea_vta_overlap(vta, atlas, side)
% Calculate the overlap between (binary) VTA and atlas.
%
% Outputs:
%     vox_overlap      : number of voxels of overlap between VTA and atlas
%     mm_overlap       : volume of overlap between VTA and atlas
%     normVTAOverlap   : normalized overlap in respect of the VTA
%     normAtlasOverlap : normalized overlap in respect of the atlas
%     mm_vta          : volume of the VTA
%     mm_atlas        : volume of the atlas

% Split left/right side of the atlas
if ~exist('side', 'var')
    side = 'both';
elseif ischar(side)
    side = lower(side);
elseif isnumeric(side)
    if side==1
        side = 'right';
    elseif  side==2
        side = 'left';
    else
        side = 'both';
    end
end

tempdir = ea_getleadtempdir;

% Load VTA image
vtanii = ea_load_nii(vta);

% Make sure the vta image is really binary
threshold_vta = max(vtanii.img(:)) * 0.5;
vtanii.img = double(vtanii.img>threshold_vta);

% Calculate volume of the VTA
mm_vta = sum(vtanii.img(:)>0) * prod(ea_detvoxsize(vta));

% Write temp VTA image
vtanii.fname = [tempdir, ea_generate_uuid, '.nii'];
ea_write_nii(vtanii);

% Threshold atlas
atlasnii = ea_load_nii(atlas);

% Only calculate for one side, suppose RAS orientation
if ~strcmp(side, 'both')
    nonzeroInd = find(atlasnii.img(:));
    [xvox, yvox, zvox] = ind2sub(size(atlasnii.img), nonzeroInd);
    xyz = ea_vox2mm([xvox, yvox, zvox], atlas);
    switch side
        case {'right', 'r'}
            atlasnii.img(nonzeroInd(xyz(:,1)<0))=0; % Set left side to 0
        case {'left', 'l'}
            atlasnii.img(nonzeroInd(xyz(:,1)>0))=0; % Set right side to 0
    end
end

% Threshold atlas
threshold_atlas = max(atlasnii.img(:)) * 0.5;
atlasnii.img = double(atlasnii.img > threshold_atlas);

% Calculate volume of the atlas
mm_atlas = sum(atlasnii.img(:)>0) * prod(ea_detvoxsize(atlas));

% Write temp atlas image
atlasnii.fname = [tempdir, ea_generate_uuid, '.nii'];
ea_write_nii(atlasnii);

% Calculate overlap
[vox_overlap, normVTAOverlap, ~, normAtlasOverlap] = ea_nii_overlap(vtanii.fname, atlasnii.fname, 1, [threshold_vta, threshold_atlas]);
mm_overlap = vox_overlap * prod(ea_detvoxsize(vta));

% Cleanup
ea_delete({vtanii.fname, atlasnii.fname});
