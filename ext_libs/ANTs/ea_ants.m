function ea_ants(fixedimage, movingimage, outputname)
% Wrapper for ANTs nonlinear registration

basedir = [fileparts(mfilename('fullpath')), filesep];

if ispc
    HEADER = [basedir, 'PrintHeader.exe'];
    ANTS = [basedir, 'antsRegistration.exe'];
elseif isunix
    HEADER = [basedir, 'PrintHeader.', computer];
    ANTS = [basedir, 'antsRegistration.', computer];
end

if fileparts(movingimage)
    volumedir = [fileparts(movingimage), filesep];
else
    volumedir =['.', filesep];
end

[~, imgsize] = system([HEADER, fixedimage, '']);
imgsize = cellfun(@(x) str2double(x),strsplit(imgsize,'x'));

if any(imgsize>256)
    rigidconvergence='[1000x500x250x0,1e-6,10]';
    rigidshrinkfactors='12x8x4x2';
    rigidsoomthingssigmas='4x3x2x1vox';

    affineconvergence='[1000x500x250x0,1e-6,10]';
    affineshrinkfactors='12x8x4x2';
    affinesoomthingssigmas='4x3x2x1vox';
else
    rigidconvergence='[1000x500x250x0,1e-6,10]';
    rigidshrinkfactors='8x4x2x1';
    rigidsoomthingssigmas='3x2x1x0vox';

    affineconvergence='[1000x500x250x0,1e-6,10]';
    affineshrinkfactors='8x4x2x1';
    affinesoomthingssigmas='3x2x1x0vox';
end

rigidstage = [' --initial-moving-transform [', fixedimage, ',', movingimage, ',1]' ...
              ' --transform Rigid[0.1]' ...
              ' --metric MI[', fixedimage, ',', movingimage, ',1,32,Regular,0.25]' ...
              ' --convergence ', rigidconvergence, ...
              ' --shrink-factors ', rigidshrinkfactors, ...
              ' --smoothing-sigmas ', rigidsoomthingssigmas];

affinestage = [' --transform Affine[0.1]'...
               ' --metric MI[', fixedimage, ',', movingimage, ',1,32,Regular,0.25]' ...
               ' --convergence ', affineconvergence, ...
               ' --shrink-factors ', affineshrinkfactors ...
               ' --smoothing-sigmas ', affinesoomthingssigmas];

ea_libs_helper
system([ANTS, ' --verbose 1' ...
              ' --dimensionality 3 --float 1' ...
              ' --output [',outputname(1:end-4), ',', outputname, ']' ...
              ' --interpolation Linear' ...
              ' --use-histogram-matching 1' ...
              ' --winsorize-image-intensities [0.005,0.995]', ...
              rigidstage, affinestage]);

movefile([volumedir, outputname(1:end-4), '0GenericAffine.mat'], [volumedir, 'ct2anat.mat']);
