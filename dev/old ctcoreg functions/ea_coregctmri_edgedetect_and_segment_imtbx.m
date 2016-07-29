function varargout=ea_coregctmri_edgedetect_and_segment_imtbx(options)
% __________________________________________________________________________________
% Copyright (C) 2014 Charite University Medicine Berlin, Movement Disorders Unit
% Andreas Horn

if ischar(options) % return name of method.
    varargout{1}='Edgedetection+Segment, use ML imagetoolbox';
    if exist('edge.m','file') % check for imtbx.
    varargout{2}={'SPM8','SPM12'};
    else
        varargout{2}={};
    end
    varargout{3}=['0.9']; % suggestion for alpha-parameter.
    return
end

if ~exist([options.root,options.patientname,filesep,'c',options.prefs.prenii_unnormalized],'file')
    ea_nc_segment(options);
end

%% run standard coregistration
options.usediffmr_coregct=['c',options.prefs.prenii_unnormalized];
ea_coregctmri_edgedetect_imtbx(options);
