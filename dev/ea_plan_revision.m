function ea_plan_revision(directory,sourceels,revs)
if exist([directory,'ea_reconstruction_backup.mat'],'file')
    copyfile( [directory,'ea_reconstruction_backup.mat'],[directory,'ea_reconstruction.mat']);
end
load([directory,'ea_reconstruction.mat']);
save([directory,'ea_reconstruction_backup.mat'],'reco');
for rev=1:length(revs)
    options=ea_getptopts(directory);

    load([directory,'ea_reconstruction.mat']);
    %    if ~isfield(reco,'acpc')
    %        options.prefs.reco.saveACPC=1;
    %        [coords_mm,trajectory,markers,elmodel,manually_corrected,coords_acpc]=ea_load_reconstruction(options);
    %        options.hybridsave=1;
    %        ea_save_reconstruction(coords_mm,trajectory,markers,elmodel,manually_corrected,options);
    %        load([directory,'ea_reconstruction.mat']);
    %    end

    copymarkers=reco.native.markers(sourceels(rev));
    % head
    cfg.acmcpc=2;
    cfg.mapmethod=0;

    cfg.xmm=[copymarkers.head(1)]; cfg.ymm=[copymarkers.head(2)]; cfg.zmm=[copymarkers.head(3)];
    fid=ea_native2acpc(cfg,{directory});
    newmarkers.head=fid.WarpedPointACPC+revs{rev}(1,:);
    % after transform back to MNI
    cfg.xmm=[newmarkers.head(1)]; cfg.ymm=[newmarkers.head(2)]; cfg.zmm=[newmarkers.head(3)];
    fid=ea_acpc2mni(cfg,{directory});
    newmarkers.head=fid.WarpedPointMNI;

    % tail
    cfg.xmm=[copymarkers.tail(1)]; cfg.ymm=[copymarkers.tail(2)]; cfg.zmm=[copymarkers.tail(3)];
    fid=ea_native2acpc(cfg,{directory});
    newmarkers.tail=fid.WarpedPointACPC+revs{rev}(2,:);

    % after transform back to MNI
    cfg.xmm=[newmarkers.tail(1)]; cfg.ymm=[newmarkers.tail(2)]; cfg.zmm=[newmarkers.tail(3)];
    fid=ea_acpc2mni(cfg,{directory});
    newmarkers.tail=fid.WarpedPointMNI;

    reco.mni.markers(end+1).head=newmarkers.head;
    reco.mni.markers(end).tail=newmarkers.tail;
    normtrajvector=(reco.mni.markers(end).tail-reco.mni.markers(end).head)./norm(reco.mni.markers(end).tail-reco.mni.markers(end).head);
    orth=null(normtrajvector)*(options.elspec.lead_diameter/2);
    reco.mni.markers(end).x=reco.mni.markers(end).head+orth(:,1)';
    reco.mni.markers(end).y=reco.mni.markers(end).head+orth(:,2)'; % corresponding points in reality

    reco=rmfield(reco,'native');
    save([directory,'ea_reconstruction.mat'],'reco');
    options.hybridsave=1;
    options.sides=1:length(reco.mni.markers);
    [reco.mni.coords_mm,reco.mni.trajectory,reco.mni.markers]=ea_resolvecoords(reco.mni.markers,options,0);
    delete([directory,'ea_reconstruction.mat']);
    ea_save_reconstruction(reco.mni.coords_mm, reco.mni.trajectory, reco.mni.markers, reco.props(1).elmodel, 1, options)
end

