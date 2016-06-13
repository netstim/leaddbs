function ea_prepare_hardi(options)
    % build HARDI in new way:
    
    if ~exist([options.root,options.patientname,filesep,options.prefs.HARDI],'file');
        hr=ea_DPS_nifti_to_hardi([options.root,options.patientname,filesep,options.prefs.dti],...
            [options.root,options.patientname,filesep,options.prefs.bvec],...
            [options.root,options.patientname,filesep,options.prefs.bval]);
        
        mrstruct_write(hr,[options.root,options.patientname,filesep,options.prefs.HARDI])
    end