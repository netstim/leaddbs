function [XYZ_mm, XYZ_src_vx] = ea_map_coords(varargin)
% This version of map_coords is based on Ged Ridgway's version but is
% optimized for usage in Lead-DBS. Especially, it supports ANTs and
% interpolations in high dimensional warping with SPM.

%  % high-dimensional warping / y_ deformation field:
%  % from target voxel space to source world and voxel space
%    [XYZ_mm XYZ_src_vx] = map_coords(XYZ_vx, '', 'y_img.nii', src);
%
% Ged Ridgway (drc.spm at gmail.com)

try
    XYZ_vx=varargin{1};
end
try
    trg=varargin{2};
end
try
    xfrm=varargin{3};
end
try
    src=varargin{4};
end

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

        % check if ANTs has been used here:
        directory=[fileparts(xfrm),filesep];
        if nargin<5
            whichnormmethod=ea_whichnormmethod(directory);
        else
            whichnormmethod=varargin{5};
        end
        if ismember(whichnormmethod,ea_getantsnormfuns)
            
            [~,fn]=fileparts(xfrm);
            if ~isempty(strfind(fn,'inv'))
                useinverse=1;
            else
                useinverse=0;
            end
            V=spm_vol(trg);
            
            % voxel to mm
            XYZ_mm_beforetransform=V(1).mat*XYZ_vx;
            
            % RAS to LPS
            XYZ_mm_beforetransform(1,:)=-XYZ_mm_beforetransform(1,:);
            XYZ_mm_beforetransform(2,:)=-XYZ_mm_beforetransform(2,:);
            
            % normalization
            XYZ_mm=ea_ants_applytransforms_to_points(directory,XYZ_mm_beforetransform(1:3,:)',useinverse)';
            XYZ_mm=[XYZ_mm;ones(1,size(XYZ_mm,2))];
            
            % LPS to RAS
            XYZ_mm(1,:)=-XYZ_mm(1,:);
            XYZ_mm(2,:)=-XYZ_mm(2,:);
        elseif ismember(whichnormmethod,ea_getfslnormfuns)
            
            [~,fn]=fileparts(xfrm);
            if ~isempty(strfind(fn,'inv'))
                useinverse=1;
            else
                useinverse=0;
            end
            V=spm_vol(trg);
            XYZ_mm_beforetransform=V(1).mat*XYZ_vx;
            prefs=ea_prefs('');
            hdr=cbiReadNiftiHeader(V.fname);
            [~,preniibase]=fileparts(prefs.gprenii);

            if ~useinverse
                ea_fsl_gencoordwarpfile(1, [directory,preniibase,'Composite.nii'],[directory,preniibase,'CompositeCoords.nii'], hdr);
                XYZ_vx_trsf = fslApplyWarpCoords(XYZ_vx,ea_detvoxsize(V(1).mat),1, [directory,preniibase,'CompositeCoords.nii'], hdr);
            else
                ea_fsl_gencoordwarpfile(1, [directory,preniibase,'InverseComposite.nii'],[directory,preniibase,'InverseCompositeCoords.nii'], hdr);
                XYZ_vx_trsf = fslApplyWarpCoords(XYZ_vx,ea_detvoxsize(V(1).mat),1, [directory,preniibase,'CompositeCoords.nii'], hdr);
            end
        else
            XYZ_mm = hdw_trgvx2srcmm(XYZ_vx, xfrm);
        end
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

XYZ_mm=XYZ_mm(1:3,:);
try
    XYZ_src_vx=XYZ_src_vx(1:3,:);
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

Vx=spm_vol([hdwim,',1,1']);
Vy=spm_vol([hdwim,',1,2']);
Vz=spm_vol([hdwim,',1,3']);

for i = 1:size(coord, 2)
    ind = inds(:, i);
    coord(1,i)=spm_sample_vol(Vx,ind(1),ind(2),ind(3),1);
    coord(2,i)=spm_sample_vol(Vy,ind(1),ind(2),ind(3),1);
    coord(3,i)=spm_sample_vol(Vz,ind(1),ind(2),ind(3),1);
end
