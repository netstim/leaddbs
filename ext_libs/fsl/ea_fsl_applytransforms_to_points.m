function outcoords = ea_fsl_applytransforms_to_points(varargin)
% Non-linear point mapping between anat and mni using FSL(img2imgcoord)
% For linear normalization, please use ea_img2imgcoord directly
% mm to mm

directory=fullfile(varargin{1},filesep);
incoords=varargin{2};

% Inverse here means REALLY inverse: map from mni coords to anat coords
if nargin == 3
    useinverse = varargin{3};
else
    useinverse = 0;
end

if nargin == 4
    transform = varargin{4};
else
    transform = '';
end

options.prefs=ea_prefs(fileparts(directory));
[~, warpprefix] = ea_niifileparts(options.prefs.gprenii);

if useinverse
    src = [ea_getearoot, 'templates', filesep, 'mni_hires.nii'];
    dest = [directory, options.prefs.prenii_unnormalized];
    if isempty(transform)
        transform = [directory, warpprefix, 'InverseWarpFiled.nii'];
    end
else
    src = [directory, options.prefs.prenii_unnormalized];
    dest = [ea_getearoot, 'templates', filesep, 'mni_hires.nii'];
    if isempty(transform)
        transform = [directory, warpprefix, 'WarpFiled.nii'];
    end
end

outcoords = ea_img2imgcoord(incoords, src, dest, transform, 'n');
