function varargout=ea_normalize_reslicepretra(options)

% __________________________________________________________________________________
% Copyright (C) 2014 Charite University Medicine Berlin, Movement Disorders Unit
% Andreas Horn

if ischar(options) % return name of method.
    varargout{1}='Reslice Preoperative MRI (for manual fusion with CT)';
    varargout{2}={'SPM8','SPM12'};
    return
end
V=spm_vol([options.root,options.patientname,filesep,options.prefs.prenii_unnormalized]);

dim=V.mat(logical(eye(4)));
dim=abs(dim(1:3));
if any(dim>0.7)
ea_reslice_nii([options.root,options.patientname,filesep,options.prefs.prenii_unnormalized],[options.root,options.patientname,filesep,options.prefs.prenii_unnormalized],[1 1 1],1);
end
