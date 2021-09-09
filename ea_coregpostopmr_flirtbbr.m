function varargout = ea_coregpostopmr_flirtbbr(options, fixed, moving, out)
% Wrapper for FSL BBR registration of post-op MRI

if ischar(options) % return name of method.
    varargout{1} = 'FLIRT BBR (Greve and Fischl 2009)';
    return
end

ea_flirtbbr(fixed, moving, out, options.prefs.mrcoreg.writeoutcoreg);
