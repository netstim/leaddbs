function ea_coregpostopmr_brainsfit(options, fixed, moving, out)
% Wrapper for BRAINSFit registration of post-op MRI

ea_brainsfit(fixed, moving, out, options.prefs.mrcoreg.writeoutcoreg);
