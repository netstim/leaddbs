function ea_atlas2stl(atlasnames,ofn)

if ~iscell(atlasnames)
    ea_error('Please specify atlas(es) in a cellstring');
end
cnt=1;
earoot=ea_getearoot;
for atl=1:length(atlasnames)
    
    load([earoot,'atlases',filesep,atlasnames{atl},filesep,'atlas_index.mat']);
    for mesh=1:numel(atlases.fv)
        cfv(cnt)=atlases.fv{mesh};
        cnt=cnt+1;
    end
end
cfv=ea_concatfv(cfv);

ea_stlwrite(ofn,cfv);