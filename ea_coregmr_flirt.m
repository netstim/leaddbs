function ea_coregmr_flirt(options, fixed, moving, out)
% Wrapper for FLIRT FLIRT registration

ea_flirt(fixed, moving, out, options.prefs.mrcoreg.writeoutcoreg);
