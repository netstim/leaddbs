function ea_pat2ls(uipatdir,handles)




    
atlasset=get(handles.atlassetpopup,'String');
atlasset=atlasset{get(handles.atlassetpopup,'Value')};


cfv(1)=ea_atlas2ply({atlasset},[uipatdir,filesep,'export',filesep,'ply',filesep,'anatomy.ply']);
try % this is DBS specific.
    cfv(2)=ea_electrode2ply([uipatdir,filesep],1,handles);
    cfv(3)=ea_electrode2ply([uipatdir,filesep],2,handles);
    cfv=ea_concatfv(cfv);
    ea_patch2ply([uipatdir,filesep,'export',filesep,'ply',filesep,'combined_scene.ply'],cfv.vertices',cfv.faces',cfv.facevertexcdata');
end


