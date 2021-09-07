function ea_coregmr_ants(options,fixed, moving, out, refine)
% uses ANTs to coregister MRIs.
if ~exist('refine','var')
    refine=0;
end

ea_ants(fixed, moving, out, options.prefs.mrcoreg.writeoutcoreg, {}, refine, options);
