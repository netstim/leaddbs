function cfv=ea_atlas2ply(atlasnames,ofn)

if ~iscell(atlasnames)
    ea_error('Please specify atlas(es) in a cellstring');
end
cnt=1;
earoot=ea_getearoot;
for atl=1:length(atlasnames)

    load([ea_space([],'atlases'),atlasnames{atl},filesep,'atlas_index.mat']);
    for side=1:2
    for mesh=1:length(atlases.names)
        cfv(cnt).vertices=atlases.fv{mesh,side}.vertices;
        cfv(cnt).faces=atlases.fv{mesh,side}.faces;
        cfv(cnt).facevertexcdata=atlases.cdat{mesh,side}';
        if isempty(cfv(cnt).facevertexcdata) % fiber atlas
            cfv(cnt).facevertexcdata=repmat(atlases.colors(mesh),size(cfv(cnt).vertices,1),1);
        end
        cnt=cnt+1;
    end
    end
end

cfv=ea_concatfv(cfv);


cfv=ea_mapcolvert2ind(cfv);
cfv.faces=[cfv.faces(:,2),cfv.faces(:,1),cfv.faces(:,3)];
ea_patch2ply(ofn,cfv.vertices',cfv.faces',cfv.facevertexcdata');
