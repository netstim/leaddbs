function ea_coregpostopmr_flirt(options, fixed, moving, out)
% Wrapper for FSL FLIRT registration of post-op MRI

ea_flirt(fixed, moving, out, options.prefs.mrcoreg.writeoutcoreg);
