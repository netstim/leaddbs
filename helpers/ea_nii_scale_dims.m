function nii_scale_dims(fnms, scale)
%Change image resolution: useful for faster classroom demos or artificially interpolated reconstructions
% fnms  : file name[s] of image[s] (optional)
% scale : size scaling factor: 2 doubles resolution, 0.5 halves resolution, "0.5 0.5 1" halves in plane
%License
% Created by Chris Rorden Sept 2014
% This function uses code from k-wave http://www.k-wave.org/license.php
% It inherits the GNU Lesser General Public License
%Examples
% nii_scale_dims; %use GUI
% nii_scale_dims('ToF.nii',2); %double resolution
% nii_scale_dims('ToF.nii',[0.5 0.5 1]); %half in plane resolution
% nii_scale_dims(strvcat('a.nii','b.nii')),0.25;

if ~exist('fnms','var')
	fnms = spm_select(inf,'image','Select image[s] for NaN removal');
end
if ~exist('scale','var')
	answer = inputdlg('Scaling factor ("0.5"=half size, "2 2 1"=double in plane)' ,'Settings',1,{'2'});
	scale=str2num(answer{1}); %#ok<ST2NM>
end
if ~exist('imresize3','file') && (min(scale(:)) < 1)
    fprintf('warning: imresize3 not found, downsampled image may have aliasing\n');
    fprintf('     https://blogs.mathworks.com/steve/2017/01/16/aliasing-and-image-resizing-part-3/');
end
if numel(scale) < 3
	if numel(scale) < 2
		scale(2) = scale(1);
	end;
	scale(3) = scale(2);
end
for i=1:size(fnms,1)
    fnm = deblank(fnms(i,:));
    fnm = unGzSub(fnm);
    hdr = spm_vol(fnm);
    %if (mod(hdr(1).dim(1),2) ~= 0) || (mod(hdr(1).dim(2),2) ~= 0)
    %    error('%s requires images to have an even number of rows and columns (not %d %d)\n',mfilename, hdr(1).dim(1), hdr(1).dim(2));
    %end
    img = spm_read_vols(hdr);
    hdrOut = hdr(1);
    hdrOut.dim(1) = round(hdr(1).dim(1)*scale(1));
    hdrOut.dim(2) = round(hdr(1).dim(2)*scale(2));
    hdrOut.dim(3) = round(hdr(1).dim(3)*scale(3));

    hdrOut.mat(1:3, 1:3) = hdr(1).mat(1:3, 1:3)*[1/scale(1) 0 0; 0 1/scale(2) 0; 0 0 1/scale(3)];
    % http://stackoverflow.com/questions/12520152/resizing-3d-matrix-image-in-matlab
    % monotonic, use the methods '*linear', '*cubic', or '*nearest'.
    [pth nm ext] = spm_fileparts(fnm);
    hdrOut.fname = fullfile(pth, ['z' nm ext]);
    fprintf('Resizing %d volumes from %dx%dx%d to %dx%dx%d\n', size(img,4), hdr(1).dim(1), hdr(1).dim(2), hdr(1).dim(3), ...
        hdrOut.dim(1), hdrOut.dim(2), hdrOut.dim(3) );
    for vol=1:size(img,4)
        hdrOut.n(1)=vol;
        if exist('imresize3','file')
            %see https://nbviewer.jupyter.org/urls/dl.dropbox.com/s/s0nw827nc4kcnaa/Aliasing.ipynb
            %downsample with anti-aliasing filter  https://en.wikipedia.org/wiki/Image_scaling
            if min(scale) < 1
                imgOut = imresize3(img,hdrOut.dim(1:3),'Antialiasing',true);
            else
                imgOut = imresize3(img,hdrOut.dim(1:3));
            end
        else
            imgOut  = resizeSub(img(:, :, :, vol), hdrOut.dim(1:3),'*cubic');
        end;
        spm_write_vol(hdrOut, imgOut);
    end;
end
%end nii_scale_dims() - local functions follow

function fnm = unGzSub(fnm);
[~,~,ext] = spm_fileparts(fnm);
if ~strcmpi(ext,'.gz'), return; end;
fnmGz = fnm;
fnm = char(gunzip(fnmGz));
delete(fnmGz); %fsl does not allow extensions .nii .nii.gz for same image
%end unGzSub()

function mat_rs = resizeSub(varargin)
%RESIZE     Resize a matrix.
% DESCRIPTION:
%       Resize a matrix to a given size using interp2 (2D) or interp3
%       (3D).
%       Use interpolation to redivide the [0,1] interval into Nx, Ny, Nz
%       voxels, where 0 is the center of first voxel, and 1 is the center
%       of the last one.
% USAGE:
%       mat_rs = resize(mat, new_size)
%       mat_rs = resize(mat, new_size, interp_mode)
%
% INPUTS:
%       mat         - matrix to resize
%       new_size    - desired matrix size in elements given by [Nx, Ny] in
%                     2D and [Nx, Ny, Nz] in 3D. Here Nx is the number of
%                     elements in the row direction, Ny is the number of
%                     elements in the column direction, and Nz is the
%                     number of elements in the depth direction.
% OPTIONAL INPUTS:
%       interp_mode - interpolation mode used by interp2 and interp3
%                     (default = '*linear')
% OUTPUTS:
%       mat_rs      - resized matrix
%disp('Resizing matrix...');
% assign the matrix input
mat = varargin{1};
% check for interpolation mode input
if nargin == 2
    interp_mode = '*linear';
elseif nargin ~= 3
    error('incorrect number of inputs');
else
    interp_mode = varargin{3};
end
% check inputs
if numDimSub(mat) ~= length(varargin{2})
    error('resolution input must have the same number of elements as data dimensions');
end
switch numDimSub(mat)
case 2
    % extract the original number of pixels from the size of the matrix
    [Nx_input, Ny_input] = size(mat);
    % extract the desired number of pixels
    Nx_output = varargin{2}(1);
    Ny_output = varargin{2}(2);
    % update command line status
    %disp(['  input grid size: ' num2str(Nx_input) ' by ' num2str(Ny_input) ' elements']);
    %disp(['  output grid size: ' num2str(Nx_output) ' by ' num2str(Ny_output) ' elements']);
    % check the size is different to the input size
    if Nx_input ~= Nx_output || Ny_input ~= Ny_output
        % resize the input matrix to the desired number of pixels
        mat_rs = interp2(0:1/(Ny_input - 1):1, (0:1/(Nx_input - 1):1)', mat, 0:1/(Ny_output - 1):1, (0:1/(Nx_output - 1):1)', interp_mode);
    else
        mat_rs = mat;
    end
case 3
    % extract the original number of pixels from the size of the matrix
    [Nx_input, Ny_input, Nz_input] = size(mat);
    % extract the desired number of pixels
    Nx_output = varargin{2}(1);
    Ny_output = varargin{2}(2);
    Nz_output = varargin{2}(3);
    % update command line status
    %disp(['  input grid size: ' num2str(Nx_input) ' by ' num2str(Ny_input) ' by ' num2str(Nz_input) ' elements']);
    %disp(['  output grid size: ' num2str(Nx_output) ' by ' num2str(Ny_output) ' by ' num2str(Nz_output) ' elements']);
    % create normalised plaid grids of current discretisation
    [x_mat, y_mat, z_mat] = ndgrid((0:Nx_input-1)/(Nx_input-1), (0:Ny_input-1)/(Ny_input-1), (0:Nz_input-1)/(Nz_input-1));
    % create plaid grids of desired discretisation
    [x_mat_interp, y_mat_interp, z_mat_interp] = ndgrid((0:Nx_output-1)/(Nx_output-1), (0:Ny_output-1)/(Ny_output-1), (0:Nz_output-1)/(Nz_output-1));
    % compute interpolation; for a matrix indexed as [M, N, P], the
    % axis variables must be given in the order N, M, P
    mat_rs = interp3(y_mat, x_mat, z_mat, mat, y_mat_interp, x_mat_interp, z_mat_interp, interp_mode);
otherwise
    error('input matrix must be 2 or 3 dimensional');
end
%end numDimSub()

function N = numDimSub(A)
N = ndims(A);
%end numDimSub()
