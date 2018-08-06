function outcoords = ea_fsl_apply_normalization_to_points(varargin)
% Non-linear point mapping between anat and mni using FSL(img2imgcoord)
% For linear transformation, please use ea_fsl_img2imgcoord directly
% mm to mm

directory=fullfile(varargin{1},filesep);
incoords=varargin{2};

% INVERSE means REALLY inverse: map from mni coords to anat coords
if nargin == 3
    inversemap = varargin{3};
else
    inversemap = 0;
end

if nargin == 4
    transform = varargin{4};
else
    transform = '';
end

options.prefs=ea_prefs(fileparts(directory));
[~, warpprefix] = ea_niifileparts(options.prefs.gprenii);

if inversemap
    src = [ea_space,'t2.nii'];
    dest = [directory, options.prefs.prenii_unnormalized];
    if isempty(transform)
        transform = [directory, warpprefix, 'WarpField.nii'];
    end
else
    src = [directory, options.prefs.prenii_unnormalized];
    dest = [ea_space,'t2.nii'];
    if isempty(transform)
        transform = [directory, warpprefix, 'InverseWarpField.nii'];
    end
end

outcoords = ea_fsl_img2imgcoord(incoords, src, dest, transform, 'n');
