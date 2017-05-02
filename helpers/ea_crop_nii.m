function ea_crop_nii(varargin)
% This function crops bounding-box of niftis to non-zeros or non-nan values.
% usage: ea_crop_nii(filename, [prefix (default: overwrite), threshstring -
% can be 'nn' for non-Nan or 'nz' for non-zero (default), cleannan, interp]).
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

if isempty(prefix)
    prefix = 'tmp';
end

if nargin > 2 % exclude nans/zeros
    nstring = varargin{3};
else
    nstring = 'nz';
end

if nargin > 3
    cleannan = varargin{4};
else
    cleannan = 0;
end

if nargin > 4
    interp = varargin{5};
else
    interp = 0;
end

V = spm_vol(filename);

[bb,vox] = ea_spm_get_bbox(V, nstring);

bb = increasebb(bb);

if any(vox<0) || any(V.mat(logical(eye(4)))<0) %Seems to be true for FSL but not other norm. methods: vox = -0.7 0.7 0.7
    ea_reslice_nii(filename,filename,abs(vox),0,[],2,[],[],0); % last zero is to not use SPM. this has shown to not work when voxel sizes are negative
    V = spm_vol(filename);
    [bb, vox] = ea_spm_get_bbox(V, nstring);
    bb = increasebb(bb);
end

dist = diff(bb); % check for weird zero bbs in small files.

if all(dist)
    
    % create sn structure:
    sn.VG = V;
    sn.VF = V;
    sn.Affine = eye(4);
    sn.Tr = [];
    sn.flags = [];
    
    % create ropts structure:
    ropts.preserve = 0;
    ropts.bb = bb;
    ropts.vox = vox;
    ropts.interp = interp;
    ropts.wrap = [0,0,0];
    ropts.prefix = prefix;
    try
        spm_write_sn([filename,',1'], sn,ropts);
    catch % init SPM first
        spm_jobman('initcfg');
        spm_write_sn([filename,',1'], sn,ropts);
    end
end

[pth, fname, ext] = fileparts(filename);
output = fullfile(pth,[prefix, fname, ext]);
if strcmp(prefix, 'tmp')
    movefile(output, filename);
    output = filename;
end

if cleannan
    nii = ea_load_nii(output);
    nii.img(abs(nii.img)<0.01) = nan; % reduce noise in originally zero compartments.
    ea_write_nii(nii);
end

if wasgz
    gzip(output);
    delete(output);
    if ~strcmp(prefix, 'tmp')  % delete unzipped file in non-overwrite mode
        delete(filename);
    end
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
