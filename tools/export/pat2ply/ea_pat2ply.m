function ea_pat2ply(uipatdir,handles,target)

if exist('target','var')
    viewsets=load([ea_getearoot,'helpers',filesep,'export',filesep,'ea_exportviews']);
    atlasset=viewsets.(target).atlas;
else
    target=[];
    atlasset=get(handles.atlassetpopup,'String');
    atlasset=atlasset{get(handles.atlassetpopup,'Value')};
end

exportFolder = [uipatdir,filesep,'export',filesep,'ply'];
ea_mkdir(exportFolder);

if isfolder(fullfile(ea_space, 'atlases', atlasset))
    cfv(1)=ea_atlas2ply({atlasset},[exportFolder,filesep,'anatomy.ply'],target);
    cnt = 1;
else
    cnt = 0;
end

options=ea_detsides(ea_getptopts(uipatdir));
for side=options.sides
    cfv(1+cnt)=ea_electrode2ply(uipatdir,side,handles);
    cnt=cnt+1;
end

cfvel=ea_concatfv(cfv(2:end));
plywrite([exportFolder,filesep,'combined_electrodes.ply'],cfvel.faces,cfvel.vertices,cfvel.facevertexcdata,repmat(100,size(cfvel.facevertexcdata,1),1));

cfv=ea_concatfv(cfv);
plywrite([exportFolder,filesep,'combined_scene.ply'],cfv.faces,cfv.vertices,cfv.facevertexcdata,repmat(100,size(cfv.facevertexcdata,1),1));
