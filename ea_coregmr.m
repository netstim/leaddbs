function [finas]=ea_coregmr(options,automan)

if options.modality==2 % in CT imaging, coregistration is done elsewhere.
    return
end

whichnormmethod=ea_whichnormmethod([options.root,options.patientname,filesep]);
if strcmp(whichnormmethod,'ea_normalize_ants')
    doreslice=1;
else
    doreslice=0;
end



switch options.coregmr.method
    case 1 % SPM
        ea_coregmr_spm(options,automan,doreslice);
        return
    case 2 % ANTs
        ea_coregmr_ants(options);
        return
    case 3 % Brainsfit
        ea_coregmr_brainsfit(options);
        return
    case 4 % Hybrid SPM & ANTs
        ea_coregmr_spm(options,automan,0); % dont use doreslice here to refrain for doing two interpolations.
        ea_coregmr_ants(options);
    case 5 % Hybrid SPM & Brainsfit
        ea_coregmr_spm(options,automan,0); % dont use doreslice here to refrain for doing two interpolations.
        ea_coregmr_brainsfit(options);
    case 6 % Do nothing
        return
end

