function ea_openpatdir(handles)

outfolder=get(handles.patdir_choosebox,'String');

if ismember(outfolder,{'No Patient Selected','Choose Patient Directory'})
    msgbox('Please set the working directory first!', 'Error','error');
    return;
end

if strcmp(outfolder(1:8),'Multiple')
   msgbox('Multiple patients selected');
   return;
end

if ismac
    system(['open ', outfolder]);
elseif isunix
    system(['xdg-open ', outfolder]);
elseif ispc
    system(['explorer ', outfolder]);
end

cd(outfolder);