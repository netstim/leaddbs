function ea_save_reconstruction(coords_mm,trajectory,markers,elmodel,manually_corrected,options)

reco.props.elmodel=elmodel;
reco.props.manually_corrected=manually_corrected;

if options.native
    reco.native.coords_mm=coords_mm;
    reco.native.trajectory=trajectory;
    reco.native.markers=markers;
    save([options.root,options.patientname,filesep,'ea_reconstruction'],'reco');
    if isfield(options,'hybridsave');
        ea_dispt('Warping fiducials to template space');
        
        ea_reconstruction2mni(options);
        if options.prefs.reco.saveACPC
            ea_dispt('Mapping fiducials to AC/PC space');
            ea_reconstruction2acpc(options);
        end
        ea_checkswap_lr(options); % PaCER support, right could be left and vice versa.
    end
else
    reco.mni.coords_mm=coords_mm;
    reco.mni.trajectory=trajectory;
    reco.mni.markers=markers;
    save([options.root,options.patientname,filesep,'ea_reconstruction'],'reco');
    
    if isfield(options,'hybridsave');
        try
            ea_dispt('Warping fiducials to native space');
            ea_reconstruction2native(options);
            if options.prefs.reco.saveACPC
                ea_dispt('Mapping fiducials to AC/PC space');
                ea_reconstruction2acpc(options);
            end
        end
    end
end


function ea_checkswap_lr(options)
options.native=0; % this can only be done in MNI space.
[coords_mm,trajectory,markers,elmodel,manually_corrected]=ea_load_reconstruction(options);

if mean(coords_mm{1}(:,1))<mean(coords_mm{2}(:,1)) % RL swapped
    % swap RL:
    options.hybridsave=1;
    ncoords_mm{1}=coords_mm{2};    ncoords_mm{2}=coords_mm{1};
    ntrajectory{1}=trajectory{2};    ntrajectory{2}=trajectory{1};
    nmarkers(1)=markers(2); nmarkers(2)=markers(1);
    
    ea_save_reconstruction(ncoords_mm,ntrajectory,nmarkers,elmodel,manually_corrected,options);
end

