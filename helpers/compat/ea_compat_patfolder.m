function ea_compat_patfolder(options)

directory=[options.root,options.patientname,filesep];

% move anatomical images
try
    movefile([directory,'anat.nii'],[directory,'anat_t2.nii']);
end

% move T2 reference image
try
    movefile([ea_getearoot,'templates',filesep,'mni_hires.nii'],[ea_space,'t2.nii']);
end

% move legacy ANTs warps
if ismember(ea_whichnormmethod(directory),ea_getantsnormfuns)
    ea_conv_antswarps(directory);
end

% move normalized Fibertracts
ea_conv_wftr(options);
