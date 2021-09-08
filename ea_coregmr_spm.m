function ea_coregmr_spm(options, fixed, moving, out, doreslice)
% Wrapper for SPM registration of post-op MRI

costfuns={'nmi','mi','ecc','ncc'};
cfundo=[2,1];
for costfun=1:length(cfundo)
    if costfun==length(cfundo) % only at final stage apply refining if set
        ea_spm_coreg(options, moving, fixed, ...
            costfuns{cfundo(costfun)},doreslice,{''},options.prefs.mrcoreg.writeoutcoreg);
    else % dont reslice, dont refine (not last pass).
        ea_spm_coreg(options, moving, fixed, ...
            costfuns{cfundo(costfun)},0,{''},options.prefs.mrcoreg.writeoutcoreg);
    end
    disp(['*** Coregistration pass (',costfuns{cfundo(costfun)},') completed.']);
end

if doreslice
    try movefile(strrep(moving, [filesep 'anat' filesep], [filesep 'anat' filesep 'r']), out); end
end

