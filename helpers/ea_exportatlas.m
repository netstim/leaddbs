function ea_exportatlas(~,~,cmd,handles)

% get atlas name
atlname=get(handles.atlassetpopup,'String');
atlname=atlname(get(handles.atlassetpopup,'Value'));

[fn,pth]=uiputfile(['*.',lower(cmd)],'Save atlas export as...',[atlname{1},'.',lower(cmd)]);

switch lower(cmd)
    case 'stl'
        ea_atlas2stl(atlname,fullfile(pth,fn));     
    case 'ply'
        ea_atlas2ply(atlname,fullfile(pth,fn));
end
