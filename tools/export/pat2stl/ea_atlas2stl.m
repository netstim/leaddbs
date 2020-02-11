function cfv=ea_atlas2stl(atlasnames,ofn,target)

if ~exist('target','var')
    target=[];
end

if ~iscell(atlasnames)
    ea_error('Please specify atlas(es) in a cellstring');
end


cnt=1;
for atl=1:length(atlasnames)
    load([ea_space([],'atlases'),atlasnames{atl},filesep,'atlas_index.mat']);
    presets=atlases.presets(1).show;
    if ~isempty(target)
        viewsets=load([ea_getearoot,'helpers',filesep,'export',filesep,'ea_exportviews']);
        views=viewsets.(target).plyview;
        presets=resolveviews(views(1).structures,atlases);
    end
    for side=1:2
        for mesh=presets
            cfv(cnt).vertices=atlases.fv{mesh,side}.vertices;
            cfv(cnt).faces=atlases.fv{mesh,side}.faces;
            cfv(cnt).facevertexcdata=repmat(atlases.colors(mesh),size(cfv(cnt).faces,1),1);
            if isempty(cfv(cnt).facevertexcdata) % fiber atlas
                cfv(cnt).facevertexcdata=repmat(atlases.colors(mesh),size(cfv(cnt).faces,1),1);
            end
            cnt=cnt+1;
        end
    end
end

cfv=ea_concatfv(cfv);
cfv=ea_mapcolvert2face(cfv);

ea_stlwrite(ofn,cfv,'FACECOLOR',cfv.facevertexcdata);
