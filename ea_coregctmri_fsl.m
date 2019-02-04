function varargout=ea_coregctmri_fsl(options)
% This function uses FSL to register postop-CT to preop-MR.
% __________________________________________________________________________________
% Copyright (C) 2019 Charite University Medicine Berlin, Movement Disorders Unit
% Andreas Horn

if ischar(options) % return name of method.
    varargout{1}='FSL FLIRT';
    varargout{2}={'SPM8','SPM12'};
    varargout{3}=['nan']; % suggestion for alpha-parameter.
    return
end

disp('Coregistering postop CT to preop MRI...');
ea_flirt([options.root,options.patientname,filesep,options.prefs.prenii_unnormalized],...
          [options.root,options.patientname,filesep,options.prefs.rawctnii_unnormalized],...
          [options.root,options.patientname,filesep,options.prefs.ctnii_coregistered]);
disp('Coregistration done.');

