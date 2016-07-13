function ea_pat2stl(uipatdir,handles)


    
atlasset=get(handles.atlassetpopup,'String');
atlasset=atlasset{get(handles.atlassetpopup,'Value')};


cfv(1)=ea_atlas2stl({atlasset},[uipatdir,filesep,'export',filesep,'stl',filesep,'anatomy.stl']);
%try % this is DBS specific.
    cfv(2)=ea_electrode2stl([uipatdir,filesep],1,handles);
    cfv(3)=ea_electrode2stl([uipatdir,filesep],2,handles);
    cfv=ea_concatfv(cfv);
    ea_stlwrite([uipatdir,filesep,'export',filesep,'stl',filesep,'combined_scene.stl'],cfv,'FACECOLOR',cfv.facevertexcdata);
%end


