function cfv=ea_atlas2stl(atlasnames,ofn,target,outputsinglefile)

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

if outputsinglefile
    cfv=ea_concatfv(cfv);
    cfv=ea_mapcolvert2face(cfv);
    ea_stlwrite(ofn,cfv,'FACECOLOR',cfv.facevertexcdata);
else
    for i=1:length(presets)
        pth = fileparts(ofn);
        if isempty(pth)
            pth = '.';
        end

        fname = [pth, filesep, atlases.labels{1}{presets(i)}, '_Right.stl'];
        cfv(i)=ea_mapcolvert2face(cfv(i));
        ea_stlwrite(fname,cfv(i),'FACECOLOR',cfv(i).facevertexcdata);

        fname = [pth, filesep, atlases.labels{1}{presets(i)}, '_Left.stl'];
        cfv(i+length(presets))=ea_mapcolvert2face(cfv(i+length(presets)));
        ea_stlwrite(fname,cfv(i+length(presets)),'FACECOLOR',cfv(i+length(presets)).facevertexcdata);
    end
end


function presets=resolveviews(structures,atlases)
atlasnames=ea_stripext(atlases.names);
presets=find(ismember(atlasnames,structures));
