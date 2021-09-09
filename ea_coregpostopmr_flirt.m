function varargout = ea_coregpostopmr_flirt(options, fixed, moving, out)
% Wrapper for FSL FLIRT registration of post-op MRI

if ischar(options) % return name of method.
    varargout{1} = 'FLIRT (Jenkinson 2001 & 2002)';
    return
end

ea_flirt(fixed, moving, out, options.prefs.mrcoreg.writeoutcoreg);
