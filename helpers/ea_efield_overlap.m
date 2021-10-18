function [overlap, normOverlap, efieldSum]  = ea_efield_overlap(vta, atlas, side)
% Calculate the overlap between VTA efield and atlas.
%
% Outputs:
%     overlap          : overlap between VTA efield and atlas (sum of the efield values within the overlap region)
%     normOverlap      : normalized overlap in respect of the VTA efield
%     efieldSum        : sum of values within the VTA efield

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

% Threshold the efield to avoid leak overlap
prefs = ea_prefs;
efieldthreshold = prefs.machine.vatsettings.horn_ethresh*10^3;
vtanii.img(vtanii.img<=efieldthreshold) = 0;
efieldSum = sum(vtanii.img(:));

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

% Write temp atlas image
atlasnii.fname = [tempdir, ea_generate_uuid, '.nii'];
ea_write_nii(atlasnii);

% Reslice atlas using VTA as reference
ea_fsl_reslice(atlasnii.fname, vta, atlasnii.fname, 'trilinear', 0);

% Calculate overlap
atlasnii = ea_load_nii(atlasnii.fname);
threshold = max(atlasnii.img(:)) * 0.5;
atlasnii.img = double(atlasnii.img > threshold);
overlap = sum(vtanii.img(:) .* atlasnii.img(:));
normOverlap = overlap/efieldSum;
if isnan(normOverlap)
    normOverlap = 0;
end

% Cleanup
ea_delete(atlasnii.fname);
