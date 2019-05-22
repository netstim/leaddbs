function cfv=ea_atlas2ply(atlasnames,ofn,target)

if ~exist('target','var')
    target=[];
end
if ~iscell(atlasnames)
    ea_error('Please specify atlas(es) in a cellstring');
end
cnt=1;
earoot=ea_getearoot;

for atl=1:length(atlasnames)
    load([ea_space([],'atlases'),atlasnames{atl},filesep,'atlas_index.mat']);
    showwhat=atlases.presets(1).show;
    if ~isempty(target)
        viewsets=load([ea_getearoot,'helpers',filesep,'export',filesep,'ea_exportviews']);
        views=viewsets.(target).plyview;
        showwhat=resolveviews(views(1).structures,atlases);
    end
    load([ea_space([],'atlases'),atlasnames{atl},filesep,'atlas_index.mat']);
    for side=1:2

        for mesh=showwhat
            cfv(cnt).vertices=atlases.fv{mesh,side}.vertices;
            cfv(cnt).faces=atlases.fv{mesh,side}.faces;
            cfv(cnt).facevertexcdata=repmat(atlases.colors(mesh),1,size(cfv(cnt).vertices,1));
            if isempty(cfv(cnt).facevertexcdata) % fiber atlas
                cfv(cnt).facevertexcdata=repmat(atlases.colors(mesh),size(cfv(cnt).vertices,1),1);
            end
            cnt=cnt+1;
        end
    end
end

if any(size(cfv(1).facevertexcdata)==1) % convert from indexed to rgb colors.
    try
        jetlist=atlases.colormap;
    catch
        jetlist=parula;
    end
    for entry=1:length(cfv)
        if size(cfv(entry).facevertexcdata,1)==1
            cfv(entry).facevertexcdata=cfv(entry).facevertexcdata';
        end
        cfv(entry).facevertexcdata=jetlist(cfv(entry).facevertexcdata,:);
    end
    %fv(f).facevertexcdata=

end




cfv=ea_concatfv(cfv);

if ~isfield(atlases,'colormap')
    atlases.colormap=jet;
end

cfv=ea_mapcolvert2ind(cfv,atlases.colormap);
cfv.faces=[cfv.faces(:,2),cfv.faces(:,1),cfv.faces(:,3)];
%ea_patch2ply(ofn,cfv.vertices',cfv.faces',cfv.facevertexcdata');
[pth]=fileparts(ofn);
if ~exist(pth,'dir')
    mkdir(pth);
end
plywrite(ofn,cfv.faces,cfv.vertices,cfv.facevertexcdata)



function showwhat=resolveviews(structures,atlases)
atlasnames=ea_stripext(atlases.names);
showwhat=find(ismember(atlasnames,structures));

