function ea_coregpostopmr_ants(options,fixed, moving, out, refine)
% Wrapper for ANTs registration of post-op MRI

if ~exist('refine','var')
    refine=0;
end

ea_ants(fixed, moving, out, options.prefs.mrcoreg.writeoutcoreg, {}, refine, options);
