function ea_crop_nii(varargin)
% This function crops bounding-box of niftis to non-zeros or non-nan values.
% usage: ea_crop_nii(filename, [prefix (default: overwrite), threshstring -
% can be 'nn' for non-Nan or 'nz' for non-zero (default), interp, rim]).
% __________________________________________________________________________________
% Copyright (C) 2015 Charite University Medicine Berlin, Movement Disorders Unit
% Andreas Horn

filename = varargin{1};

if strcmp(filename(end-2:end), '.gz')
    wasgz = 1;
    gunzip(filename);
    filename = filename(1:end-3);
else
    wasgz = 0;
end

if nargin > 1
    prefix = varargin{2};
else
    prefix = '';
end

% if isempty(prefix)
%     prefix = 'tmp';
% end



if nargin > 3
    interp = varargin{4};
else
    interp = 0;
end


% load nifti
nii=ea_load_nii(filename);

if ~any(nii.img(:))
    return
end

if all(isnan(nii.img(:)))
    return
end

if nargin > 2 % exclude nans/zeros
    nstring = varargin{3};
    if isempty(nstring)
        if any(isnan(nii.img(:)))
            nstring='nn';
        else
            nstring = 'nz';
        end
    end
else
    if any(isnan(nii.img(:)))
        nstring='nn';
    else
        nstring = 'nz';
    end
end

switch nstring
    case 'nz'
        [xx,yy,zz]=ind2sub(size(nii.img),find(nii.img~=0));
    case 'nn'
        [xx,yy,zz]=ind2sub(size(nii.img),find(~isnan(nii.img)));
    otherwise
        [xx,yy,zz]=ind2sub(size(nii.img),find(nii.img~=num2str(nstring)));
end
if nargin>4
    rim=varargin{5}; % go to 3 mm, assuming isotropic image. No worries if not isotropic though, then rim will be a bit asymmetrical.
else    
    rim=round(2/mean(nii.voxsize)); % go to 3 mm, assuming isotropic image. No worries if not isotropic though, then rim will be a bit asymmetrical.
end
bbim=[min(xx),max(xx)
    min(yy),max(yy)
    min(zz),max(zz)];

offset=[min(xx)-rim-1
    min(yy)-rim-1
    min(zz)-rim-1];

try
switch nstring
    case 'nz'
        X=zeros((bbim(1,2)-bbim(1,1))+1+2*rim,(bbim(2,2)-bbim(2,1))+1+2*rim,(bbim(3,2)-bbim(3,1))+1+2*rim);
    case 'nn'
        X=nan((bbim(1,2)-bbim(1,1))+1+2*rim,(bbim(2,2)-bbim(2,1))+1+2*rim,(bbim(3,2)-bbim(3,1))+1+2*rim);
end
catch
    keyboard
end
tmat=eye(4);
tmat(1:3,4)=offset(:,1);



X((rim+1):end-rim,(rim+1):end-rim,(rim+1):end-rim)=nii.img(bbim(1,1):bbim(1,2),bbim(2,1):bbim(2,2),bbim(3,1):bbim(3,2));

nii.mat=nii.mat*tmat;
nii.img=X;
nii.dim=size(nii.img);

[pth,fn,ext]=fileparts(nii.fname);
nii.fname=fullfile(pth,[prefix,fn,ext]);
if isempty(prefix) % make sure nifti is fully overwritten with new copy from RAM.
    delete(filename);
end
ea_write_nii(nii);


if wasgz
    gzip(nii.fname);
    delete(nii.fname);
end



function [BB,vx] = ea_spm_get_bbox(V, thr, premul)
% Compute volume's bounding box, for full field of view or object bounds
% FORMAT [BB,vx] = spm_get_bbox(V, thr)
% V   - mapped image volume(s) (from spm_vol) or filename (empty for GUI)
% thr - threshold, such that BB contains voxels with intensities > thr
%       or strings 'nz', 'nn', fv', for non-zero, non-NaN, or field of view
%       where 'fv' (the default) uses only the image's header information.
%
% BB  - a [2 x 3] array of the min and max X, Y, and Z coordinates {mm},
%       i.e. BB = [minX minY minZ; maxX maxY maxZ].
% vx  - a [1 x 3] vector of voxel dimensions {mm}.
%__________________________________________________________________________
% Copyright (C) 2011-2013 Wellcome Trust Centre for Neuroimaging

% Ged Ridgway
% $Id: spm_get_bbox.m 5398 2013-04-12 12:37:00Z ged $

% Undocumented expert options:
% V           - can be a 4D @nifti object (but not 5D), image-based BBs
%               will be computed using "all" along the 4th dimension.
% thr = 'old' - reproduce spm_write_sn/bbvox_from_V (and elsewhere)
% premul      - a matrix that premultiplies V.mat, as used in spm_orthviews

%-Get an SPM volume structure
%--------------------------------------------------------------------------
if nargin < 1 || isempty(V)
    [V, sts] = spm_select(1, 'image', 'Select Image');
    if ~sts, error('Must select an image'), end
end
if ischar(V), V = spm_vol(V); end

%-Get volume structure from @nifti object if given
%--------------------------------------------------------------------------
if isa(V, 'nifti')
    V = spm_vol(V.dat.fname); % (potentially a struct array of volumes)
end

%-Compute voxel dimensions (for compatibility with bbvox_from_V)
%--------------------------------------------------------------------------
P = spm_imatrix(V(1).mat);
vx = P(7:9);
% the above agrees with sqrt(sum(V.mat(1:3,1:3).^2)) for simple rotations,
% and seems more appropriate if there are reflections and/or skews.
% Note that spm_imatrix(diag([-1 1 1 1])) is [-1 1 1] as expected.

%-Compute bounding box
%--------------------------------------------------------------------------
if nargin < 2 || isempty(thr) || strcmpi(thr, 'fv')
    % overall field-of-view bounding box from header information
    d = V(1).dim;
    corners = [
        1    1    1    1
        1    1    d(3) 1
        1    d(2) 1    1
        1    d(2) d(3) 1
        d(1) 1    1    1
        d(1) 1    d(3) 1
        d(1) d(2) 1    1
        d(1) d(2) d(3) 1
        ]';
    XYZ = V(1).mat(1:3, :) * corners;
elseif strcmpi(thr, 'old')
    % code from spm_write_sn/bbvox_from_V (and other places)
    % NB: main difference is that vx(1)<0 gives descending BB(:,1),
    % shouldn't be used if V.mat contains rotations or skews.
    o  = V(1).mat\[0 0 0 1]';
    o  = o(1:3)';
    BB = [-vx.*(o-1) ; vx.*(V(1).dim(1:3)-o)];
    if exist('premul', 'var')
        warning('spm_get_bbox:old_and_premul', 'old method ignores premul')
    end
else
    % image-based bounding box using voxel intensities
    img = spm_read_vols(V);
    if ischar(thr)
        switch lower(thr)
            case 'nn'  % non-NaN, though include +/- Inf in computation
                img = ~isnan(img);
            case 'nz'  % special case of non-zero (rather than > 0)
                img = ~isnan(img) & img ~= 0;
            otherwise
                error('Unknown threshold type %s', thr)
        end
    else
        % treat thr as numeric threshold
        img = img > thr;
    end
    if ndims(img) == 4
        img = all(img, 4);
    end
    if nnz(img) == 0
        warning('spm_get_bbox:nothing', ...
            'Threshold leaves no voxels, returning full field of view');
        if exist('premul', 'var')
            [BB,vx] = spm_get_bbox(V, 'fv', premul);
        else
            [BB,vx] = spm_get_bbox(V, 'fv');
        end
        return
    else
        img = find(img); % (clears img to save memory)
        [X Y Z] = ind2sub(V(1).dim, img);
        XYZ = V(1).mat(1:3, :) * [X Y Z ones(size(X))]';
    end
end

if ~exist('BB', 'var') % exists already if 'old' case chosen above
    if exist('premul', 'var')
        XYZ = premul(1:3, :) * [XYZ; ones(1, size(XYZ, 2))];
    end
    BB = [
        min(XYZ, [], 2)'
        max(XYZ, [], 2)'
        ];
end


function bb=increasebb(bb)

for dim=1:3
    dbb=bb(:,dim);
    
    [~,minpos]=min(dbb);
    [~,maxpos]=max(dbb);
    dbb(minpos)=dbb(minpos)-5;
    dbb(maxpos)=dbb(maxpos)+5;
    
    bb(:,dim)=dbb;
end
