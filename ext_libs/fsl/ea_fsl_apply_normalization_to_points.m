function outcoords = ea_fsl_apply_normalization_to_points(varargin)
% Non-linear point mapping between anat and mni using FSL(img2imgcoord)
% For linear transformation, please use ea_fsl_img2imgcoord directly
% mm to mm

options = ea_getptopts(varargin{1}); % varargin{1} is patient folder
incoords = varargin{2};

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

if inversemap
    src = [ea_space, 't2.nii'];
    dest = options.subj.preopAnat.(options.subj.AnchorModality).coreg;
    if isempty(transform)
        transform = [options.subj.norm.transform.forwardBaseName, 'fnirt.nii'];
    end
else
    src = options.subj.preopAnat.(options.subj.AnchorModality).coreg;
    dest = [ea_space, 't2.nii'];
    if isempty(transform)
        transform = [options.subj.norm.transform.inverseBaseName, 'fnirt.nii'];
    end
end

outcoords = ea_fsl_img2imgcoord(incoords, src, dest, transform, 'nonlinear');
