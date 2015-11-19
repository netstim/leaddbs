function ea_brainsfit(fixedVolume, movingVolume, outputVolume)
% Wrapper for BRAINSFit

if fileparts(fixedVolume)
    volumedir = [fileparts(fixedVolume), filesep];
else
    volumedir =['.', filesep];
end

fixparams = [' --useRigid --useAffine' ...
             ' --samplingPercentage 0.005' ...
             ' --removeIntensityOutliers 0.005' ...
             ' --interpolationMode Linear' ...
             ' --outputTransform ', volumedir, 'ct2anat.txt'];

if exist([volumedir, 'ct2anat.txt'],'file')
    fixparams = [fixparams, ...
                 ' --initializeTransformMode Off' ...
                 ' --initialTransform ', volumedir, 'ct2anat.txt'];
else
	fixparams = [fixparams, ' --initializeTransformMode useGeometryAlign']; 
end

basename = [fileparts(mfilename('fullpath')), filesep, 'BRAINSFit'];

if ispc
    BRAINSFit = [basename, '.exe '];
elseif isunix
    BRAINSFit = [basename, '.', computer, ' '];
end

ea_libs_helper
system([BRAINSFit, fixparams, ...
        ' --fixedVolume ' , fixedVolume, ...
        ' --movingVolume ', movingVolume, ...
        ' --outputVolume ', outputVolume]);
delete([volumedir, 'ct2anat_Inverse.h5']);