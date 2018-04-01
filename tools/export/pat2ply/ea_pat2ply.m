function ea_pat2ply(uipatdir,handles)




    
atlasset=get(handles.atlassetpopup,'String');
atlasset=atlasset{get(handles.atlassetpopup,'Value')};


cfv(1)=ea_atlas2ply({atlasset},[uipatdir,filesep,'export',filesep,'ply',filesep,'anatomy.ply']);
try % this is DBS specific.
    cfv(2)=ea_electrode2ply([uipatdir,filesep],1,handles);
    cfv(3)=ea_electrode2ply([uipatdir,filesep],2,handles);
    cfv=ea_concatfv(cfv);
    plywrite([uipatdir,filesep,'export',filesep,'ply',filesep,'combined_scene.ply'],cfv.faces,cfv.vertices,cfv.facevertexcdata,repmat(100,size(cfv.facevertexcdata,1),1));
    %write_ply(cfv.vertices',cfv.faces',[uipatdir,filesep,'export',filesep,'ply',filesep,'combined_scene.ply']);
    %ea_patch2ply([uipatdir,filesep,'export',filesep,'ply',filesep,'combined_scene.ply'],cfv.vertices',cfv.faces',cfv.facevertexcdata');
end


