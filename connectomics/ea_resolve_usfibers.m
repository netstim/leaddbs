function fibers=ea_resolve_usfibers(options,fibers)

if options.lc.struc.ft.upsample.factor>1 && options.lc.struc.ft.upsample.how==0
    Vb0upsampled=ea_open_vol([options.root,options.patientname,filesep,options.prefs.b0]);
    % get fibers to mm
    tempfib=[fibers(:,1:3)';ones(1,size(fibers,1))]; % vox in upsampled space
    tempfib=Vb0upsampled.mat*tempfib; % mm
    
    % pop stashed (uninterpolated) files:
    % pop dti / b0 files:
    movefile([options.root,options.patientname,filesep,ea_stripext(options.prefs.dti)],...
        [options.root,options.patientname,filesep,options.prefs.dti]);
    movefile([options.root,options.patientname,filesep,ea_stripext(options.prefs.b0)],...
        [options.root,options.patientname,filesep,options.prefs.b0]);
    Vb0regular=ea_open_vol([options.root,options.patientname,filesep,options.prefs.b0]);
    tempfib=Vb0regular.mat\tempfib; % vox in regular space 
    fibers(:,1:3)=tempfib(1:3,:)';
end