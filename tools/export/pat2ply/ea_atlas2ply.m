function cfv=ea_atlas2ply(atlasnames,ofn,target,outputsinglefile)

if ~iscell(atlasnames)
    ea_error('Please specify atlas(es) in a cellstring');
end

if ~exist('target','var')
    target=[];
end

if ~exist('outputsinglefile','var')
    outputsinglefile = 1;
end

cnt = 1;
for atl=1:length(atlasnames)
    load([ea_space([],'atlases'),atlasnames{atl},filesep,'atlas_index.mat']);
    if ~isempty(target)
        viewsets=load([ea_getearoot,'helpers',filesep,'export',filesep,'ea_exportviews']);
        views=viewsets.(target).plyview;
        presets=resolveviews(views(1).structures,atlases);
    else
        if isfield(atlases,'presets')
            presets=atlases.presets(1).show;
        else
            presets=1:length(atlases.names);
        end
    end
    for side=1:2
        for i=1:length(presets)
            cfv(cnt).vertices=atlases.roi{presets(i),side}.fv.vertices;
            cfv(cnt).faces=atlases.roi{presets(i),side}.fv.faces;
            cfv(cnt).facevertexcdata=repmat(atlases.roi{presets(i),side}.color,size(cfv(cnt).vertices,1),1);
            cnt = cnt + 1;
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
end

if ~isfield(atlases,'colormap')
    atlases.colormap=jet;
end

if outputsinglefile
    cfv=ea_concatfv(cfv);
    cfv=ea_mapcolvert2ind(cfv,atlases.colormap);
    cfv.faces=[cfv.faces(:,2),cfv.faces(:,1),cfv.faces(:,3)];
    [pth]=fileparts(ofn);
    if ~exist(pth,'dir')
        mkdir(pth);
    end
    plywrite(ofn,cfv.faces,cfv.vertices,cfv.facevertexcdata);
else
    for i=1:length(presets)
        pth = fileparts(ofn);
        if isempty(pth)
            pth = '.';
        end

        fname = [pth, filesep, atlases.labels{1}{presets(i)}, '_Right.ply'];
        cfv(i)=ea_mapcolvert2ind(cfv(i),atlases.colormap);
        cfv(i).faces=[cfv(i).faces(:,2),cfv(i).faces(:,1),cfv(i).faces(:,3)];
        plywrite(fname,cfv(i).faces,cfv(i).vertices,cfv(i).facevertexcdata);

        fname = [pth, filesep, atlases.labels{1}{presets(i)}, '_Left.ply'];
        cfv(i+length(presets))=ea_mapcolvert2ind(cfv(i+length(presets)),atlases.colormap);
        cfv(i+length(presets)).faces=[cfv(i+length(presets)).faces(:,2),cfv(i+length(presets)).faces(:,1),cfv(i+length(presets)).faces(:,3)];
        plywrite(fname,cfv(i+length(presets)).faces,cfv(i+length(presets)).vertices,cfv(i+length(presets)).facevertexcdata);
    end
end


function presets=resolveviews(structures,atlases)
atlasnames=ea_stripext(atlases.names);
presets=find(ismember(atlasnames,structures));
