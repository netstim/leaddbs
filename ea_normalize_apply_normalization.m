function varargout=ea_normalize_apply_normalization(options)
% This is a function that simply applies normalization parameters (y_ea_normparams.nii) to the
% unnormalized files and coregisters postop MRI to preop MRI if
% MRI-modality is used.
%
% __________________________________________________________________________________
% Copyright (C) 2014 Charite University Medicine Berlin, Movement Disorders Unit
% Andreas Horn


if ischar(options) % return name of method.
    varargout{1}='(Re-)apply (priorly) estimated normalization.';
    varargout{2}={'SPM8','SPM12'};
    return
end

disp('(Re)-applying priorly estimated normalizations...');
if options.modality==1 % MR
ea_coreg(options,'auto');
end
ea_apply_normalization(options);


function [flirtmat spmvoxmat fslvoxmat] = worldmat2flirtmat(worldmat, src, trg)
%worldmat2flirtmat: convert NIfTI world (mm) coordinates matrix to flirt
%
% Example:
%  [flirtmat spmvoxmat fslvoxmat] = worldmat2flirtdmat(worldmat, src, trg);
%
% See also: flirtmat2worldmat, flirtmat_write

% Copyright 2009 Ged Ridgway <ged.ridgway gmail.com>

if ischar(src)
    src = nifti(src);
end
if ischar(trg)
    trg = nifti(trg);
end

spmvoxmat = inv(src.mat) * worldmat * trg.mat;
addone = eye(4); addone(:, 4) = 1;
fslvoxmat = inv(addone) * spmvoxmat * addone;
trgscl = nifti2scl(trg);
srcscl = nifti2scl(src);
flirtmat = inv( srcscl * fslvoxmat * inv(trgscl) );



function flirtmat_write(fname, mat)
%flirtmat_write: save a 4-by-4 matrix to a file as handled by flirt -init
% Example:
%  flirtmat_write(fname, affinemat)
% See also: flirtmat_read, flirtmat2worldmat, worldmat2flirtmat

% Copyright 2009 Ged Ridgway <ged.ridgway gmail.com>

fid = fopen(fname, 'w');
% Note that transpose is needed to go from MATLAB's column-major matrix to
% writing the file in a row-major way (see also flirtmat_read)
str = sprintf('%f  %f  %f  %f\n', mat');
fprintf(fid, str(1:end-1)); % drop final newline
fclose(fid);




function [worldmat spmvoxmat fslvoxmat] = flirtmat2worldmat(flirtmat, src, trg)
%flirtmat2worldmat: convert saved flirt matrix to NIfTI world coords matrix
% flirt matrix is from text file specified in "flirt -omat mat.txt" command
% world matrix maps from NIfTI world coordinates in target to source. Note:
% mat.txt contains a mapping from source to target in FSL's *scaled* coords
% which are not NIfTI world coordinates, and note src-trg directionality!
% worldmat from this script reproduces "img2imgcoord -mm ...".
%
% The script can also return a matrix to map from target to source voxels
% in MATLAB/SPM's one-based convention, or in FSL's zero-based convention
%
% Example:
%  [worldmat spmvoxmat fslvoxmat] = flirtmat2worldmat(flirtmat, src, trg);
%
% See also: worldmat2flirtmat, flirtmat_read, flirtmat_write

% Copyright 2009 Ged Ridgway <ged.ridgway gmail.com>

if ischar(src)
    src = nifti(src);
end
if ischar(trg)
    trg = nifti(trg);
end

if ischar(flirtmat)
    flirtmat = flirtmat_read(flirtmat);
end

% src = inv(flirtmat) * trg
% srcvox = src.mat \ inv(flirtmat) * trg.mat * trgvox
% BUT, flirt doesn't use src.mat, only absolute values of the
% scaling elements from it,
% AND, if images are not radiological, the x-axis is flipped, see:
%  https://www.jiscmail.ac.uk/cgi-bin/webadmin?A2=ind0810&L=FSL&P=185638
%  https://www.jiscmail.ac.uk/cgi-bin/webadmin?A2=ind0903&L=FSL&P=R93775
trgscl = nifti2scl(trg);
srcscl = nifti2scl(src);
fslvoxmat = inv(srcscl) * inv(flirtmat) * trgscl;

% AND, Flirt's voxels are zero-based, while SPM's are one-based...
addone = eye(4); addone(:, 4) = 1;
spmvoxmat = addone * fslvoxmat * inv(addone);

worldmat = src.mat * spmvoxmat * inv(trg.mat);


%%
function scl = nifti2scl(N)
% not sure if this is always correct with rotations in mat, but seems okay!
scl = diag([sqrt(sum(N.mat(1:3,1:3).^2)) 1]);
if det(N.mat) > 0
    % neurological, x-axis is flipped, such that [3 2 1 0] and [0 1 2 3]
    % have the same *scaled* coordinates:
    xflip = diag([-1 1 1 1]);
    xflip(1, 4) = N.dat.dim(1)-1; % reflect about centre
    scl = scl * xflip;
end

function resize_img(imnames, Voxdim, BB, ismask)
%  resize_img -- resample images to have specified voxel dims and BBox
% resize_img(imnames, voxdim, bb, ismask)
%
% Output images will be prefixed with 'r', and will have voxel dimensions
% equal to voxdim. Use NaNs to determine voxdims from transformation matrix
% of input image(s).
% If bb == nan(2,3), bounding box will include entire original image
% Origin will move appropriately. Use world_bb to compute bounding box from
% a different image.
%
% Pass ismask=true to re-round binary mask values (avoid
% growing/shrinking masks due to linear interp)
%
% See also voxdim, world_bb

% Based on John Ashburner's reorient.m
% http://www.sph.umich.edu/~nichols/JohnsGems.html#Gem7
% http://www.sph.umich.edu/~nichols/JohnsGems5.html#Gem2
% Adapted by Ged Ridgway -- email bugs to drc.spm@gmail.com

% This version doesn't check spm_flip_analyze_images -- the handedness of
% the output image and matrix should match those of the input.

% Check spm version:
if exist('spm_select','file') % should be true for spm5
    spm5 = 1;
elseif exist('spm_get','file') % should be true for spm2
    spm5 = 0;
else
    ea_error('Can''t find spm_get or spm_select; please add SPM to path')
end

spm_defaults;

% prompt for missing arguments
if ( ~exist('imnames','var') || isempty(char(imnames)) )
    if spm5
        imnames = spm_select(inf, 'image', 'Choose images to resize');
    else
        imnames = spm_get(inf, 'img', 'Choose images to resize');
    end
end
% check if inter fig already open, don't close later if so...
Fint = spm_figure('FindWin', 'Interactive'); Fnew = [];
if ( ~exist('Voxdim', 'var') || isempty(Voxdim) )
    Fnew = spm_figure('GetWin', 'Interactive');
    Voxdim = spm_input('Vox Dims (NaN for "as input")? ',...
        '+1', 'e', '[nan nan nan]', 3);
end
if ( ~exist('BB', 'var') || isempty(BB) )
    Fnew = spm_figure('GetWin', 'Interactive');
    BB = spm_input('Bound Box (NaN => original)? ',...
        '+1', 'e', '[nan nan nan; nan nan nan]', [2 3]);
end
if ~exist('ismask', 'var')
    ismask = false;
end
if isempty(ismask)
    ismask = false;
end

% reslice images one-by-one
vols = spm_vol(imnames);
for V=vols'
    % (copy to allow defaulting of NaNs differently for each volume)
    voxdim = Voxdim;
    bb = BB;
    % default voxdim to current volume's voxdim, (from mat parameters)
    if any(isnan(voxdim))
        vprm = spm_imatrix(V.mat);
        vvoxdim = vprm(7:9);
        voxdim(isnan(voxdim)) = vvoxdim(isnan(voxdim));
    end
    voxdim = voxdim(:)';

    mn = bb(1,:);
    mx = bb(2,:);
    % default BB to current volume's
    if any(isnan(bb(:)))
        vbb = world_bb(V);
        vmn = vbb(1,:);
        vmx = vbb(2,:);
        mn(isnan(mn)) = vmn(isnan(mn));
        mx(isnan(mx)) = vmx(isnan(mx));
    end

    % voxel [1 1 1] of output should map to BB mn
    % (the combination of matrices below first maps [1 1 1] to [0 0 0])
    mat = spm_matrix([mn 0 0 0 voxdim])*spm_matrix([-1 -1 -1]);
    % voxel-coords of BB mx gives number of voxels required
    % (round up if more than a tenth of a voxel over)
    imgdim = ceil(mat \ [mx 1]' - 0.1)';

    % output image
    VO            = V;
    [pth,nam,ext] = fileparts(V.fname);
    VO.fname      = fullfile(pth,['r' nam ext]);
    VO.dim(1:3)   = imgdim(1:3);
    VO.mat        = mat;
    VO = spm_create_vol(VO);
    spm_progress_bar('Init',imgdim(3),'reslicing...','planes completed');
    for i = 1:imgdim(3)
        M = inv(spm_matrix([0 0 -i])*inv(VO.mat)*V.mat);
        img = spm_slice_vol(V, M, imgdim(1:2), 1); % (linear interp)
        if ismask
            img = round(img);
        end
        spm_write_plane(VO, img, i);
        spm_progress_bar('Set', i)
    end
    spm_progress_bar('Clear');
end
% call spm_close_vol if spm2
if ~spm5
    spm_close_vol(VO);
end
if (isempty(Fint) && ~isempty(Fnew))
    % interactive figure was opened by this script, so close it again.
    close(Fnew);
end
disp('Done.')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function bb = world_bb(V)
%  world-bb -- get bounding box in world (mm) coordinates

d = V.dim(1:3);
% corners in voxel-space
c = [ 1    1    1    1
    1    1    d(3) 1
    1    d(2) 1    1
    1    d(2) d(3) 1
    d(1) 1    1    1
    d(1) 1    d(3) 1
    d(1) d(2) 1    1
    d(1) d(2) d(3) 1 ]';
% corners in world-space
tc = V.mat(1:3,1:4)*c;

% bounding box (world) min and max
mn = min(tc,[],2)';
mx = max(tc,[],2)';
bb = [mn; mx];