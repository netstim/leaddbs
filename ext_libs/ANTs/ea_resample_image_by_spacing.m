function ea_resample_image_by_spacing(inputimage, spacing, dosmooth, addvox, nninterp, outputimage, dimension)
% Resample image by spacing, wrapper for ANTs ResampleImageBySpacing
%
%	spacing: output image spacing, 1*3 vector
%   dosmooth: smoooth output image, default = 0
%   addvox: pad each dimension by addvox, default = 0
%   nn-interp: use NearestNeighbor interpolation or not, default = 0 (Linear interpolation)
%   dimension: image dimension, default=3

if ~exist('dosmooth', 'var')
    dosmooth = 0;
end

if ~exist('addvox', 'var')
    addvox = 0;
end

if ~exist('nn-interp', 'var')
    nninterp = 0;
end

if ~exist('outputimage', 'var')
    outputimage = inputimage;
end

if ~exist('dimension', 'var')
    V=ea_open_vol(inputimage);
    if V.volnum>1
        dimension = 4;
    else
        dimension = 3;
    end
end

dosmooth = num2str(dosmooth);
addvox = num2str(addvox);
nninterp = num2str(nninterp);
dimension = num2str(dimension);
switch dimension
    case '3'
        spacing = [num2str(spacing(1)), ' ', num2str(spacing(2)), ' ', num2str(spacing(3))];
    case '4'
        if length(spacing)==4
            spacing = [num2str(spacing(1)), ' ', num2str(spacing(2)), ' ', num2str(spacing(3)), ' ', num2str(spacing(4))];
        else
            spacing = [num2str(spacing(1)), ' ', num2str(spacing(2)), ' ', num2str(spacing(3)), ' 15'];
        end
end
basedir=[ea_getearoot,'ext_libs',filesep,'ANTs',filesep];
if ispc
    ResampleImageBySpacing = ea_path_helper([basedir, 'ResampleImageBySpacing.exe']);
else
    ResampleImageBySpacing = [basedir, 'ResampleImageBySpacing.', computer('arch')];
end

cmd = [ResampleImageBySpacing,' ',dimension, ...
                              ' ',ea_path_helper(inputimage), ...
                              ' ',ea_path_helper(outputimage), ...
                              ' ',spacing, ...
                              ' ',dosmooth, ...
                              ' ',addvox, ...
                              ' ',nninterp];

ea_libs_helper;

fprintf('\nResampling image spacing to [%s]: %s\n\n', spacing, inputimage);
if ~ispc
    system(['bash -c "', cmd, '"']);
else
    system(cmd);
end
fprintf('\n');

