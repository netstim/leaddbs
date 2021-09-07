function ea_coregmr_flirt_bbr(options, fixed, moving, out)
% uses FSL bbr instead of SPM to coregister MRIs.

ea_flirt_bbr(fixed, moving, out, options.prefs.mrcoreg.writeoutcoreg);
