function ea_export_depth_of_field(~,~,figurehandle)

[FileName,PathName] = uiputfile('LEAD_Scene.png','Save file name');
if FileName
    im=ea_depth_of_field(figurehandle);

    imwrite(im, [PathName,FileName], 'png');
end


