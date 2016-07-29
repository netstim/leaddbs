function varargout=ea_coregctmri_edgedetect_and_segment(options)
% __________________________________________________________________________________
% Copyright (C) 2014 Charite University Medicine Berlin, Movement Disorders Unit
% Andreas Horn

if ischar(options) % return name of method.
    varargout{1}='Edgedetection+Segment';
    varargout{2}={'SPM8','SPM12'};
    varargout{3}=['0.6,0.4']; % suggestion for alpha-parameter.
    return
end

if ~exist([options.root,options.patientname,filesep,'c',options.prefs.prenii_unnormalized],'file')
    ea_nc_segment(options);
end

%% run standard coregistration
options.usediffmr_coregct=['c',options.prefs.prenii_unnormalized];
ea_coregctmri_edgedetect(options);