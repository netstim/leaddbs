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
    hdr = ea_fslhd(inputimage);
    dimension = hdr.dim0;
end

if dimension==4 && length(spacing)==3
    % keep the original time step (pixdim[4])
    hdr = ea_fslhd(inputimage);
    spacing = [spacing, hdr.pixdim4];
end

dosmooth = num2str(dosmooth);
addvox = num2str(addvox);
nninterp = num2str(nninterp);
dimension = num2str(dimension);
spacing = sprintf('% f', spacing);

basedir=[ea_getearoot,'ext_libs',filesep,'ANTs',filesep];
ResampleImageBySpacing = ea_getExec([basedir, 'ResampleImageBySpacing'], escapePath = 1);


cmd = [ResampleImageBySpacing,' ',dimension, ...
                              ' ',ea_path_helper(inputimage), ...
                              ' ',ea_path_helper(outputimage), ...
                              ' ',spacing, ...
                              ' ',dosmooth, ...
                              ' ',addvox, ...
                              ' ',nninterp];

fprintf('\nResampling image spacing to [%s]: %s\n\n', spacing, inputimage);

ea_libs_helper;
ea_runcmd(cmd);

fprintf('\n');
