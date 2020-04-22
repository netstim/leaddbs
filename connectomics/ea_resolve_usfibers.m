function fibers = ea_resolve_usfibers(options, fibers)

if options.lc.struc.ft.upsample.how==0 % Dyrby et al. 2014
    % Change ZERO-BASED indexing to ONE-BASED indexing.
    fibers(:,1:3) = fibers(:,1:3) + 1;

    % Convert fibers in upsampled voxel space to mm
    ref = [options.root,options.patientname,filesep,options.prefs.b0];
    fibers(:,1:3) = ea_vox2mm(fibers(:,1:3), ref);

    % Pop stashed (uninterpolated) dti and b0:
    movefile([options.root,options.patientname,filesep,ea_stripext(options.prefs.dti)],...
        [options.root,options.patientname,filesep,options.prefs.dti]);
    movefile([options.root,options.patientname,filesep,ea_stripext(options.prefs.b0)],...
        [options.root,options.patientname,filesep,options.prefs.b0]);

    % Convert mm fibers into regular voxel space
    fibers(:,1:3) = ea_mm2vox(fibers(:,1:3), ref);
else % DSI-Studio internal upsampling
    fibers(:,1:3) = fibers(:,1:3)./ea_resolve_usfactor(options.lc.struc.ft.upsample);
    % Change ZERO-BASED indexing to ONE-BASED indexing.
    fibers(:,1:3) = fibers(:,1:3) + 1;
end
