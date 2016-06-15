function ea_atlas2stl(atlasfile,ofn)

load(atlasfile);
for mesh=1:numel(atlases.fv)
    cfv(mesh)=atlases.fv{mesh};
end
cfv=ea_concatfv(cfv);

ea_stlwrite(ofn,cfv);