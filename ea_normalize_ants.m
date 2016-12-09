function varargout=ea_normalize_ants(options)
% This is a function that normalizes both a copy of transversal and coronar
% images into MNI-space. The goal was to make the procedure both robust and
% automatic, but still, it must be said that normalization results should
% be taken with much care because all reconstruction results heavily depend
% on these results. Normalization of DBS-MR-images is especially
% problematic since usually, the field of view doesn't cover the whole
% brain (to reduce SAR-levels during acquisition) and since electrode
% artifacts can impair the normalization process. Therefore, normalization
% might be best archieved with other tools that have specialized on
% normalization of such image data.
%
% The procedure used here uses the ANTs Syn approach to map a patient's
% brain to MNI space directly. 
% __________________________________________________________________________________
% Copyright (C) 2015 Charite University Medicine Berlin, Movement Disorders Unit
% Andreas Horn


if ischar(options) % return name of method.
    varargout{1}='Advanced Normalization Tools (ANTs) SyN (Avants 2008)';
    varargout{2}={'SPM8','SPM12'};
    return
end


% ANTs nolinear registration
directory=[options.root,options.patientname,filesep];
ea_ants_nonlinear([options.earoot,'templates',filesep,'mni_hires',options.primarytemplate,'.nii'],...
                  [directory,options.prefs.prenii_unnormalized],...
                  [directory,options.prefs.gprenii]);

% Apply registration
ea_apply_normalization(options)
