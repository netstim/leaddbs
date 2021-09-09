function varargout = ea_coregpostopmr_brainsfit(options, fixed, moving, out)
% Wrapper for BRAINSFit registration of post-op MRI

if ischar(options) % return name of method.
    varargout{1} = 'BRAINSFit (Johnson 2007)';
    return
end

ea_brainsfit(fixed, moving, out, options.prefs.mrcoreg.writeoutcoreg);
