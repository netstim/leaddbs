function ea_ants_nonlinear(fixedimage, movingimage, outputimage)
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

synstage = [' --transform SyN[0.3]'...
    ' --metric MI[', fixedimage, ',', movingimage, ',1,32,Regular,0.25]' ...
    ' --convergence ', affineconvergence, ...
    ' --shrink-factors ', affineshrinkfactors ...
    ' --smoothing-sigmas ', affinesoomthingssigmas];

%           ANTS 3 -m PR[/images/subject/b0diffusionimage.nii.gz,/images/subject/T1image.nii.gz,1,2] -i 50x20x10 -o T1_to_diffusion_synants.nii.gz -t SyN[0.3] -r Gauss[3,0]

ea_libs_helper
[outputdir, outputname, ~] = fileparts(outputimage);
if outputdir
    outputbase = [outputdir, filesep, outputname];
else
    outputbase = ['.', filesep, outputname];
end
system([ANTS, ' --verbose 1' ...
    ' --dimensionality 3 --float 1' ...
    ' --output [',outputbase, ',', outputimage, ']' ...
    ' --interpolation Linear' ...
    ' --use-histogram-matching 1' ...
    ' --winsorize-image-intensities [0.005,0.995]', ...
    rigidstage, affinestage, synstage]);

%movefile([outputbase, '0GenericAffine.mat'], [volumedir, 'ct2anat.mat']);
