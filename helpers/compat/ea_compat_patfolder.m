function ea_compat_patfolder(options)

directory=[options.root,options.patientname,filesep];
% move legacy ANTs warps
ea_conv_antswarps(directory);
% move normalized Fibertracts
ea_conv_wftr(options);