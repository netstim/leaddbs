function varargout=ea_coregctmri_matlab_imreg(options)
% This function uses the Matlab image toolbox to register postop-CT to preop-MR.
% __________________________________________________________________________________
% Copyright (C) 2015 Charite University Medicine Berlin, Movement Disorders Unit
% Andreas Horn

if ischar(options) % return name of method.
    varargout{1}='BRAINSFit';
    varargout{2}={'SPM8','SPM12'};
    varargout{3}=['nan']; % suggestion for alpha-parameter.
    return
end


system([options.earoot,'dev/BRAINSFIT/BRAINSFit --fixedVolume ',[options.root,options.patientname,filesep,options.prefs.prenii_unnormalized],...
    ' --movingVolume ',[options.root,options.patientname,filesep,options.prefs.rawctnii_unnormalized],...
    ' --outputVolume ',[options.root,options.patientname,filesep,options.prefs.ctnii_coregistered],...
    ' --useRigid --useAffine --samplingPercentage 0.005 --removeIntensityOutliers 0.005 --initializeTransformMode useGeometryAlign --interpolationMode Linear']);


