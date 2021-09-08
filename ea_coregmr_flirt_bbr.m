function ea_coregmr_flirt_bbr(options, fixed, moving, out)
% Wrapper for FSL BBR linear registration

ea_flirt_bbr(fixed, moving, out, options.prefs.mrcoreg.writeoutcoreg);
