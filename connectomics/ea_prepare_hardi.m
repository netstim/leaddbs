function ea_prepare_hardi(options,redo)
    
% build HARDI in new way:
    
if ~exist([options.root,options.patientname,filesep,options.prefs.HARDI],'file') || redo;
    hr=ea_DPS_nifti_to_hardi([options.root,options.patientname,filesep,options.prefs.dti],...
        [options.root,options.patientname,filesep,options.prefs.bvec],...
        [options.root,options.patientname,filesep,options.prefs.bval]);
    
    mrstruct_write(hr,[options.root,options.patientname,filesep,options.prefs.HARDI]);
    % convert to DTD
    hr=mrstruct_read([options.root,options.patientname,filesep,options.prefs.HARDI]);
    hr=convertHARDI2DTD(hr,0,0);
    save([options.root,options.patientname,filesep,options.prefs.DTD],'-struct','hr');
    
else
    disp('HARDI found, no need to rebuild new one.');
end
    
