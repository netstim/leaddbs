function ea_coregmr_flirt(options, fixed, moving, out)
% uses FLIRT instead of SPM to coregister MRIs.

ea_flirt(fixed, moving, out, options.prefs.mrcoreg.writeoutcoreg);
