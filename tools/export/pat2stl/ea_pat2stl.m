function ea_pat2stl(uipatdir,handles)

atlasset=get(handles.atlassetpopup,'String');
atlasset=atlasset{get(handles.atlassetpopup,'Value')};

cfv(1)=ea_atlas2stl({atlasset},[uipatdir,filesep,'export',filesep,'stl',filesep,'anatomy.stl']);

cfv(2)=ea_electrode2stl([uipatdir,filesep],1,handles);
cfv(3)=ea_electrode2stl([uipatdir,filesep],2,handles);

ecfv = ea_concatfv(cfv(2:3));
ea_stlwrite([uipatdir,filesep,'export',filesep,'stl',filesep,'combined_electrodes.stl'],ecfv,'FACECOLOR',ecfv.facevertexcdata);

scfv=ea_concatfv(cfv);
ea_stlwrite([uipatdir,filesep,'export',filesep,'stl',filesep,'combined_scene.stl'],scfv,'FACECOLOR',scfv.facevertexcdata);
