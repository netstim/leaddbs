function ea_save_reconstruction(coords_mm,trajectory,markers,elmodel,manually_corrected,options)

reco.mni.coords_mm=coords_mm;
reco.mni.trajectory=trajectory;
reco.mni.markers=markers;
reco.props.elmodel=elmodel;
reco.props.manually_corrected=manually_corrected;


save([options.root,options.patientname,filesep,'ea_reconstruction'],'reco');
ea_reconstruction2native(options);