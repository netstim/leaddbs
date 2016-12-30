function varargout=ea_normalize_reslicepretra(options)

% __________________________________________________________________________________
% Copyright (C) 2014 Charite University Medicine Berlin, Movement Disorders Unit
% Andreas Horn

if ischar(options) % return name of method.
    varargout{1}='Reslice Preoperative MRI (for manual fusion with CT)';
    varargout{2}=1;
    return
end
directory=[options.root,options.patientname,filesep];
V=spm_vol([directory,options.prefs.prenii_unnormalized]);

dim=V.mat(logical(eye(4)));
dim=abs(dim(1:3));
if any(dim>0.7)
    ea_reslice_nii([directory,options.prefs.prenii_unnormalized],[directory,options.prefs.prenii_unnormalized],[0.7 0.7 0.7],0,[],1);
end
