function ea_pat2stl(hObj,~,handles)

uipatdir=getappdata(handles.leadfigure,'uipatdir');
if strcmp(uipatdir{1},'No Patient Selected')
    ea_error('Please select a patient first.');
    return
end
    
atlasset=get(handles.atlassetpopup,'String');
atlasset=atlasset{get(handles.atlassetpopup,'Value')};

for pt=1:length(uipatdir)
    mkdir([uipatdir{pt},filesep,'stlexport']);
    cfv(1)=ea_atlas2stl({atlasset},[uipatdir{pt},filesep,'stlexport',filesep,'anatomy.stl']);
    cfv(2)=ea_electrode2stl([uipatdir{pt},filesep],1);
    cfv(3)=ea_electrode2stl([uipatdir{pt},filesep],2);
    cfv=ea_concatfv(cfv);
    ea_stlwrite([uipatdir{pt},filesep,'stlexport',filesep,'combined_scene.stl'],cfv,'FACECOLOR',cfv.facevertexcdata);

end

