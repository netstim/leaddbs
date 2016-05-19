function ea_bias_field_correction(inputimage)
% Wrapper for ANTs N4BiasFieldCorrection

ea_libs_helper;

basedir = [fileparts(mfilename('fullpath')), filesep];

if ispc
    N4BiasFieldCorrection = [basedir, 'N4BiasFieldCorrection.exe'];
elseif isunix
    N4BiasFieldCorrection = [basedir, 'N4BiasFieldCorrection.', computer];
end

cmd=[N4BiasFieldCorrection, ...
    ' --image-dimensionality 3' ...
    ' --input-image ', inputimage, ...
    ' --output ', inputimage, ...
    ' --shrink-factor 4' ...
    ' --bspline-fitting [200]' ...
    ' --convergence [50x50x50x50,0.000001]'];
       
if ~ispc
    system(['bash -c "', cmd, '"']);
else
    system(cmd);
end
