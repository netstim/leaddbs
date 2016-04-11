function ea_checkforwardinv(options,forwardinv)

switch forwardinv
    case 'forward'
        V=spm_vol([options.root,options.patientname,filesep,'y_ea_normparams.nii']);
        Vmni=spm_vol([options.earoot,'templates',filesep,'mni_hires.nii']);
        
        if ~isequal(V.dim,Vmni.dim)
            
            ea_redo_inv([options.root,options.patientname,filesep],options,'forward');
        end
    case 'inverse'
        V=spm_vol([options.root,options.patientname,filesep,'y_ea_inv_normparams.nii']);
        Vanat=spm_vol([options.root,options.patientname,filesep,options.prefs.prenii_unnormalized]);
        if ~isequal(V.dim,Vanat.dim)
            ea_redo_inv([options.root,options.patientname,filesep],options,'inverse');
        end
end