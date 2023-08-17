function ea_pat2ply(uipatdir,handles,target)

if exist('target','var')
    viewsets=load([ea_getearoot,'helpers',filesep,'export',filesep,'ea_exportviews']);
    atlasset=viewsets.(target).atlas;
else
    target=[];
    atlasset=get(handles.atlassetpopup,'String');
    atlasset=atlasset{get(handles.atlassetpopup,'Value')};
end

cfv(1)=ea_atlas2ply({atlasset},[uipatdir,filesep,'export',filesep,'ply',filesep,'anatomy.ply'],target);

options=ea_detsides(ea_getptopts(uipatdir));
cnt=1;
for side=options.sides
    cfv(1+cnt)=ea_electrode2ply(uipatdir,side,handles);
    cnt=cnt+1;
end

cfvel=ea_concatfv(cfv(2:end));
plywrite([uipatdir,filesep,'export',filesep,'ply',filesep,'combined_electrodes.ply'],cfvel.faces,cfvel.vertices,cfvel.facevertexcdata,repmat(100,size(cfvel.facevertexcdata,1),1));

cfv=ea_concatfv(cfv);
plywrite([uipatdir,filesep,'export',filesep,'ply',filesep,'combined_scene.ply'],cfv.faces,cfv.vertices,cfv.facevertexcdata,repmat(100,size(cfv.facevertexcdata,1),1));
