function ea_save_reconstruction(coords_mm,trajectory,markers,elmodel,manually_corrected,options)

reco.props.elmodel=elmodel;
reco.props.manually_corrected=manually_corrected;

if options.native
    reco.native.coords_mm=coords_mm;
    reco.native.trajectory=trajectory;
    reco.native.markers=markers;
    save([options.root,options.patientname,filesep,'ea_reconstruction'],'reco');
    if ~isfield(reco,'mni')
        ea_reconstruction2mni(options);
    end
else
    reco.mni.coords_mm=coords_mm;
    reco.mni.trajectory=trajectory;
    reco.mni.markers=markers;
    save([options.root,options.patientname,filesep,'ea_reconstruction'],'reco');
    if ~isfield(reco,'native')
    ea_reconstruction2native(options);
    end
end




