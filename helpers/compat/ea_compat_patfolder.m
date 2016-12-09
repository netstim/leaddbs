function ea_compat_patfolder(options)

directory=[options.root,options.patientname,filesep];
% move legacy ANTs warps
ea_conv_antswarps(directory);
% move normalized Fibertracts
ea_conv_wftr(options);
try
movefile([directory,'anat.nii'],[directory,'anat_t2.nii']);
end
try
    movefile([ea_getearoot,'templates',filesep,'mni_hires.nii'],[ea_getearoot,'templates',filesep,'mni_hires_t2.nii']);
end