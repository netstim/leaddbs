function ea_coregpostopmr_flirtbbr(options, fixed, moving, out)
% Wrapper for FSL BBR registration of post-op MRI

ea_flirtbbr(fixed, moving, out, options.prefs.mrcoreg.writeoutcoreg);
