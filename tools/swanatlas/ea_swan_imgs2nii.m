%Author:Eduardo Alho
%date 19-02-2024
%root dir is the directory where the folder Masks is.
%This script will concatenate the original histological masks into a
%isotropic nifti file
%the resolution can be set in resolution
%Example AllMasks2Nii ('/Users/eduardoalho/Desktop/Atlas/Atlas_test/',[0.3 0.3 0.3])


function outfile=ea_swan_imgs2nii(indir, resolution)

%Get the directory name to name the nifti file
[~, structure_name, ~] = fileparts(indir);

% Read TIFF images
ims = dir(fullfile(indir, '*.tif'));
if isempty(ims)
    ims=dir(fullfile(indir,'*.png'));
end
ims(147) = [];  % Exclude specific image


% Initialize empty volume
volume = [];

% Read the first image to determine size
first_im = imread(fullfile(indir, ims(1).name));
im_size = size(first_im);
xfrom = 850:6000;

% Loop through TIFF images and stack them into volume
ea_dispercent(0,'Reading images');
for i = 1:length(ims)
    im = imread(fullfile(indir, ims(i).name));
    trans = repmat(im(:,:,1)~=255, 1, 1, 3);
    im(~trans) = 255;
    im = im(:, xfrom, :);

    % Append to volume
    if isempty(volume)
        volume = zeros([im_size(1), length(xfrom), length(ims)], 'like', first_im);
    end
    volume(:, :, i) = mean(im, 3);
    ea_dispercent(i/length(ims));
end
ea_dispercent(1,'end');

% Create NIfTI file
V = struct();
outfile = fullfile(indir, [structure_name, '.nii']);
V.fname = outfile;
V.dim = size(volume);
V.dt = [2, 0];
V.n = [1, 1];
V.descrip = 'SWANhisto';
V.pinfo = [1; 0; 352];
V.mat = load(fullfile (ea_getearoot,'tools','swanatlas','header_correct.mat')).V.mat;

% Create volume header
V = spm_create_vol(V);

% Write out NIfTI slices
ea_dispercent(0,'Creating Nifti file');
for cnt = 1:size(volume, 3)
    spm_write_plane(V, squeeze(volume(:,:,cnt)), cnt);
    ea_dispercent(cnt / size(volume, 3));
end
ea_dispercent(1, 'end');

% Reslice
ea_reslice_nii(V.fname, V.fname, resolution, [1], [0], [1], [], [], [3]);
