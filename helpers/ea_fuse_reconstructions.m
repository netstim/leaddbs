function ea_fuse_reconstructions(recos,outputfilename)
% simple function to fuse several reconstructions into one.
% recos is a cell pointing to different reconstruction iles

for r=1:length(recos)

    rec(r)=load(recos{r});
    if r == 1
        reco=rec(1).reco;
    else
        try reco.native.coords_mm=[reco.native.coords_mm,rec(r).reco.native.coords_mm]; end
        try reco.mni.coords_mm=[reco.native.coords_mm,rec(r).reco.mni.coords_mm]; end
        try reco.scrf.coords_mm=[reco.native.coords_mm,rec(r).reco.scrf.coords_mm]; end
        try reco.acpc.coords_mm=[reco.native.coords_mm,rec(r).reco.acpc.coords_mm]; end


        try reco.native.trajectory=[reco.native.trajectory,rec(r).reco.native.trajectory]; end
        try reco.mni.trajectory=[reco.mni.trajectory,rec(r).reco.mni.trajectory]; end
        try reco.scrf.trajectory=[reco.scrf.trajectory,rec(r).reco.scrf.trajectory]; end
        try reco.acpc.trajectory=[reco.acpc.trajectory,rec(r).reco.acpc.trajectory]; end

        try reco.native.markers(end+1:end+length(rec(r).reco.native.markers))=rec(r).reco.native.markers; end
        try reco.mni.markers(end+1:end+length(rec(r).reco.mni.markers))=rec(r).reco.mni.markers; end
        try reco.scrf.markers(end+1:end+length(rec(r).reco.scrf.markers))=rec(r).reco.scrf.markers; end       
        try reco.acpc.markers(end+1:end+length(rec(r).acpc.scrf.markers))=rec(r).acpc.scrf.markers; end
    end
end

save(outputfilename,'reco');
