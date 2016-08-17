function iscoreg=ea_seemscoregistered(options)

% function checks if all anatomical images in folder are coregistered or
% not.

iscoreg=1;
directory=[options.root,options.patientname,filesep];

[~,peerpresentfiles]=ea_assignpretra(options);
if exist(ea_niigz([directory,options.prefs.fa2anat]),'file')
    peerpresentfiles=[peerpresentfiles,{options.prefs.fa2anat}];
end

if length(peerpresentfiles)==1 % only one anatomical image present, nothing to coregister
    return
end

Vref=ea_open_vol([directory,peerpresentfiles{1}]);
for comp=2:length(peerpresentfiles)
    Vcomp=ea_open_vol([directory,peerpresentfiles{comp}]);
    if ~isequal(Vref.dim,Vcomp.dim)
        iscoreg=0;
        return
    end
    if ~isequal(Vref.mat,Vcomp.dim)
        iscoreg=0;
        return
    end
end