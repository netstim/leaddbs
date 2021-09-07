function ea_coregmr_brainsfit(options, fixed, moving, out)
% uses Brainsfit to coregister MRIs.

ea_brainsfit(fixed, moving, out, options.prefs.mrcoreg.writeoutcoreg);
