function ea_pat2stl(uipatdir,handles)

atlasset=get(handles.atlassetpopup,'String');
atlasset=atlasset{get(handles.atlassetpopup,'Value')};

cfv(1)=ea_atlas2stl({atlasset},[uipatdir,filesep,'export',filesep,'stl',filesep,'anatomy.stl']);

options=ea_detsides(ea_getptopts(uipatdir));
cnt=1;
for side=options.sides
    cfv(1+cnt)=ea_electrode2stl(uipatdir,side,handles);
    cnt=cnt+1;
end

ecfv = ea_concatfv(cfv(2:end));
ea_stlwrite([uipatdir,filesep,'export',filesep,'stl',filesep,'combined_electrodes.stl'],ecfv,'FACECOLOR',ecfv.facevertexcdata);

scfv=ea_concatfv(cfv);
ea_stlwrite([uipatdir,filesep,'export',filesep,'stl',filesep,'combined_scene.stl'],scfv,'FACECOLOR',scfv.facevertexcdata);
