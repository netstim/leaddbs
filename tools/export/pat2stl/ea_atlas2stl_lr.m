function cfv=ea_atlas2stl_lr(atlasnames,ofn)

if ~iscell(atlasnames)
    ea_error('Please specify atlas(es) in a cellstring');
end
earoot=ea_getearoot;
for atl=1:length(atlasnames)
    load([ea_space([],'atlases'),atlasnames{atl},filesep,'atlas_index.mat']);
    for side=1:2
        cnt=1;
        clear cfv
        for mesh=1:length(atlases.names)
            cfv(cnt).vertices=atlases.roi{mesh,side}.fv.vertices;
            cfv(cnt).faces=atlases.roi{mesh,side}.fv.faces;
            cfv(cnt).facevertexcdata=atlases.roi{mesh,side}.cdat';
            if isempty(cfv(cnt).facevertexcdata) % fiber atlas
                cfv(cnt).facevertexcdata=repmat(atlases.roi(mesh,side).color,size(cfv(cnt).faces,1),1);
            end

            cnt=cnt+1;
        end

        cfv=ea_concatfv(cfv);

        cfv=ea_mapcolvert2face(cfv);

        [pth,fn,ext]=fileparts(ofn);
        sofn=fullfile(pth,[fn,'_',num2str(side),ext]);
        ea_stlwrite(sofn,cfv,'FACECOLOR',cfv.facevertexcdata);
    end
end

