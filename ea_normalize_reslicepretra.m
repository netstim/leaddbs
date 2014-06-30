function varargout=ea_normalize_reslicepretra(options)

% __________________________________________________________________________________
% Copyright (C) 2014 Charite University Medicine Berlin, Movement Disorders Unit
% Andreas Horn

if ischar(options) % return name of method.
    varargout{1}='Reslice Preoperative MRI (for manual fusion with CT)';
    return
end

reslice_nii([options.root,options.patientname,filesep,options.prefs.prenii_unnormalized],[options.root,options.patientname,filesep,options.prefs.prenii_unnormalized],[0.5 0.5 0.5]);

