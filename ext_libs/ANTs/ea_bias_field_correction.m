function ea_bias_field_correction(inputimage)
% Wrapper for ANTs N4BiasFieldCorrection

ea_libs_helper;

basedir = [fileparts(mfilename('fullpath')), filesep];

N4BiasFieldCorrection = ea_getExec([basedir, 'N4BiasFieldCorrection'], escapePath = 1);


cmd=[N4BiasFieldCorrection, ...
    ' --image-dimensionality 3' ...
    ' --input-image ', ea_path_helper(inputimage), ...
    ' --output ', ea_path_helper(inputimage), ...
    ' --shrink-factor 4' ...
    ' --bspline-fitting [200]' ...
    ' --convergence [50x50x50x50,0.000001]'];

fprintf('\nBias field correction...\n\n')

ea_runcmd(cmd);
