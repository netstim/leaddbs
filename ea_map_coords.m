function [XYZ_mm, XYZ_src_vx] = ea_map_coords(XYZ_vx, trg, xfrm, src)
% map_coords -- map between coordinate systems/spaces
% from (target) image voxel space to world/mm space, optionally to
% registered (source) world coords, and optionally to source voxel space.
% (special case: can also map from world space to voxel space for target,
% see first example below)
%
% Registered coordinates (in world/mm space) can be due to:
%  - affine registration (e.g. reorient, realign, coreg)
%  - spatial normalisation using DCT
%    (e.g. SPM2 "Normalise" sn.mat or SPM5 "Segment" seg_sn.mat)
%  - high dimensional warping (using HDW toolbox)
%
% Coords can be passed in for n points as a 3-by-n matrix with rows for
% x, y and z components; or as a 4-by-n matrix of homogeneous coordinates
% (i.e. with a fourth row of all ones). The same will be returned.
%
% Pass empty arguments ('') to get GUI selection prompts.
% If an sn.mat is passed, trg will be ignored, and the target from the sn
% will be used; if src is unspecified, the source from the sn will be used.
% If an HDW deformation field is selected trg is ignored.
%
% Examples:
%
%  % map from world to voxel coords (special case)
%    [XYZ_mm XYZ_vx] = map_coords(XYZ_mm, img); % (XYZ_mm unaltered)
%
%  % affine:
%   % map from voxel to world coords:
%    [XYZ_mm] = map_coords(XYZ_vx, img);
%   % map to world coordinates, and then to voxel coords in source image:
%   % (note, use zero (0) to specify no transformation; empty ('') prompts)
%    [XYZ_mm XYZ_src_vx] = map_coords(XYZ_vx, trg, 0, src);
%
%  % normalisation:
%   % map from template voxel coords to DCT normalised source world coords
%    [XYZ_mm] = map_coords(XYZ_vx, '', 'blah_sn.mat');
%   % same, and then to source voxel coords as well
%    [XYZ_mm XYZ_src_vx] = map_coords(XYZ_vx, '', 'sn.mat', src);
%   % same, using source from sn.mat
%    [XYZ_mm XYZ_src_vx] = map_coords(XYZ_vx, '', 'sn.mat');
%   % to map from source to template, invert the sn.mat (e.g. using the
%   % deformations utility and then use the resulting y_ deformation field)
%
%  % unified segmentation:
%   % from template voxel space to native subject world and voxel space:
%    [XYZ_mm XYZ_src_vx] = map_coords(XYZ_vx, '', 'seg_sn.mat');
%   % from native subject voxel space to template world and voxel space:
%    [XYZ_mm XYZ_src_vx] = map_coords(XYZ_vx, '', 'seg_inv_sn.mat');
%
%  % high-dimensional warping / y_ deformation field:
%   % from target voxel space to source world and voxel space
%    [XYZ_mm XYZ_src_vx] = map_coords(XYZ_vx, '', 'y_img.nii', src);
%
% Ged Ridgway (drc.spm at gmail.com)

if nargin < 2
    error('map_coords:usage',...
        'Must specify at least coords and trg; empty [] for GUI prompt')
end

% Input coordinates
if isempty(XYZ_vx)
    XYZ_vx = spm_input('voxel coords? ', '+1', 'r', '1 1 1', [Inf 3])';
end
n = size(XYZ_vx, 2); % number of points
if size(XYZ_vx, 1) == 3
    % note to return 3*n later
    homog = false;
    % make homogeneous
    XYZ_vx = [XYZ_vx; ones(1, n)];
elseif size(XYZ_vx, 1) == 4
    homog = true;
else
    error('map_coords:dims',...
        'coord array must have 3 or 4 rows: [x;y;z] or [x;y;z;1]')
end

% Coordinate mapping
if ~exist('xfrm', 'var') || ~ischar(xfrm)
    % affine only
    if isempty(trg)
        trg = spm_select(1, 'image', 'Choose target image');
    end
    trg = spm_vol(trg);
    XYZ_mm = trg.mat * XYZ_vx;
    xfrm = []; % (now set to empty)
elseif isempty(xfrm)
    xfrm = spm_select([0 1], 'any',...
        'Select transformation file (e.g. sn.mat or HDW y_ field',...
        '', pwd, '\.(mat|img|nii)$'); % (empty if user chooses nothing)
end
if ~isempty(xfrm)
    if ~isempty(regexp(xfrm, 'sn\.mat$', 'once'))
        % DCT sn structure
        XYZ_mm = sn_trgvx2srcmm(XYZ_vx, xfrm);
    elseif ~isempty(regexp(xfrm, 'y_.*(nii|img)$', 'once'))
        % HDW transformation field
        XYZ_mm = hdw_trgvx2srcmm(XYZ_vx, xfrm);
    else
        error('map_coords:xfrm', 'unrecognised transformation file')
    end
elseif ~exist('XYZ_mm', 'var')
    error('map_coords:nothingdoing',...
        'failed to select target image or transformation, nothing to do!')
end

% Optional mapping from source world space (mm) to source voxel space
if nargout > 1
    if exist('src', 'var')
        if isempty(src)
            src = spm_select(1, 'image', 'Choose source image');
        end
        src = spm_vol(src);
        XYZ_src_vx = src.mat \ XYZ_mm;
    elseif ischar(xfrm) && ~isempty(regexp(xfrm, 'sn\.mat$', 'once'))
        load(xfrm, 'VF');
        XYZ_src_vx = VF.mat \ XYZ_mm;
    elseif exist('trg', 'var')
        % assume input actually world coords, and desired vox coord output
        XYZ_mm = XYZ_vx;
        XYZ_src_vx = trg.mat \ XYZ_mm;
    else
        error('map_coords:src_vx',...
            'source image (or sn.mat) not specified, but src_vx requested')
    end
end

% Use consistent style of coordinates between input and output
if ~homog % drop homogeneous ones from output
    XYZ_mm(4, :) = [];
    if exist('XYZ_src_vx', 'var')
        XYZ_src_vx(4, :) = [];
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function coord = sn_trgvx2srcmm(coord, matname)
% src_mm = sn_trgvx2srcmm(trg_vx, matname)
% Based on John Ashburner's get_orig_coord5.m
sn = load(matname); Tr = sn.Tr;
if numel(Tr) ~= 0 % DCT warp: trg_vox displacement
    d = sn.VG(1).dim(1:3); % (since VG may be 3-vector of TPM volumes)
    dTr = size(Tr);
    basX = spm_dctmtx(d(1), dTr(1), coord(1,:)-1);
    basY = spm_dctmtx(d(2), dTr(2), coord(2,:)-1);
    basZ = spm_dctmtx(d(3), dTr(3), coord(3,:)-1);
    for i = 1:size(coord, 2)
        bx = basX(i, :);
        by = basY(i, :);
        bz = basZ(i, :);
        tx = reshape(...
            reshape(Tr(:,:,:,1),dTr(1)*dTr(2),dTr(3))*bz',dTr(1),dTr(2) );
        ty = reshape(...
            reshape(Tr(:,:,:,2),dTr(1)*dTr(2),dTr(3))*bz',dTr(1),dTr(2) );
        tz =  reshape(...
            reshape(Tr(:,:,:,3),dTr(1)*dTr(2),dTr(3))*bz',dTr(1),dTr(2) );
        coord(1:3,i) = coord(1:3,i) + [bx*tx*by' ; bx*ty*by' ; bx*tz*by'];
    end
end
% Affine: trg_vox (possibly displaced by above DCT) to src_vox
coord = sn.VF.mat * sn.Affine * coord;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function coord = hdw_trgvx2srcmm(coord, hdwim)
% src_mm = hdw_trgvx2srcmm(trg_vx, y_hdw)
inds = coord(1:3, :);
% if abs(coord(1:3, :) - inds) > 0.01
%     warning('map_coords:hdw_rounding',...
%         'target voxel coords rounded for HDW')
% end
W = nifti(hdwim);

for i = 1:size(coord, 2)
    ind = inds(:, i);
%                % linearly interpolate: 
%     samp(1,1,1,:)=squeeze(W.dat(floor(ind(1)), floor(ind(2)), floor(ind(3)), 1, 1:3));
%     samp(1,1,2,:)=squeeze(W.dat(floor(ind(1)), floor(ind(2)), ceil(ind(3)), 1, 1:3));
%     samp(1,2,1,:)=squeeze(W.dat(floor(ind(1)), floor(ind(2)), ceil(ind(3)), 1, 1:3));
%     samp(2,1,1,:)=squeeze(W.dat(ceil(ind(1)), floor(ind(2)), floor(ind(3)), 1, 1:3));
%     samp(1,2,2,:)=squeeze(W.dat(floor(ind(1)), ceil(ind(2)), ceil(ind(3)), 1, 1:3));
%     samp(2,1,2,:)=squeeze(W.dat(ceil(ind(1)), floor(ind(2)), ceil(ind(3)), 1, 1:3));
%     samp(2,2,1,:)=squeeze(W.dat(ceil(ind(1)), ceil(ind(2)), floor(ind(3)), 1, 1:3));
%     samp(2,2,2,:)=squeeze(W.dat(ceil(ind(1)), ceil(ind(2)), ceil(ind(3)), 1, 1:3));
%     
%     rests=ind-floor(ind);
% 
%    
                coord(1:3, i) = squeeze(W.dat(ind(1), ind(2), ind(3), 1, 1:3));
end